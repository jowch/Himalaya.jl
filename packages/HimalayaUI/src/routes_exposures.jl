using HTTP, JSON3, DBInterface, Tables, Oxygen

const PRP_DETAILS = Dict{String,Dict{Symbol,Any}}()
const PRP_DETAILS_LOCK = ReentrantLock()

function prp_details(path; parse = parse_prp)
    empty = Dict{Symbol,Any}(:beam_energy_ev => nothing, :pipe_length_m => nothing, :detector => nothing)
    path isa AbstractString || return empty
    hit = lock(() -> get(PRP_DETAILS, path, nothing), PRP_DETAILS_LOCK)
    hit === nothing || return hit
    isfile(path) || return empty
    prp = try
        parse(path)
    catch err
        @warn "prp_details: could not read" path exception = err
        return empty
    end
    details = Dict{Symbol,Any}(
        :beam_energy_ev => coalesce(prp.beam_energy_ev, nothing),
        :pipe_length_m  => coalesce(prp.pipe_length_m, nothing),
        :detector       => coalesce(prp.detector, nothing),
    )
    lock(() -> PRP_DETAILS[path] = details, PRP_DETAILS_LOCK)
    details
end

function register_exposures_routes!()
    @get "/api/samples/{id}/exposures" function(req::HTTP.Request, id::Int)
        db = current_db()
        sql = "SELECT * FROM exposures WHERE sample_id = ? ORDER BY id"
        exs = sqlrows(DBInterface.execute(db, sql, [id]))
        out = map(exs) do e
            srcs = sqlrows(DBInterface.execute(db,
                "SELECT source_exposure_id, role FROM exposure_sources
                 WHERE averaged_exposure_id = ?", [Int(e.id)]))
            d = row_to_json(e)
            delete!(d, :q_step)
            delete!(d, :reviewed_at)
            delete!(d, :peak_model_version)
            d[:sources]       = rows_to_json(srcs)
            d[:image_version] = image_version_token(e.image_path)
            merge!(d, prp_details(e.prp_path))
            d
        end
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(out))
    end

    @get "/api/exposures/{id}/image" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = sqlrows(DBInterface.execute(db,
            "SELECT image_path FROM exposures WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "image not found")))

        ip = rows[1].image_path
        (ip === nothing || ip isa Missing) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "no TIFF for this image")))

        # The DB stores image_path at ingest, but the source TIFF may have been
        # moved or deleted since. Guard before load_and_lognormalize, which
        # would otherwise throw an unhandled ArgumentError → HTTP 500. Fail
        # gracefully with a 404 like the trace route does for a missing .dat
        # (issue #233, finding BE-1).
        path = String(ip)
        isfile(path) || return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "image source file not found: $path")))

        params   = HTTP.queryparams(req)
        is_thumb = get(params, "thumb", "0") == "1"

        # Compute the cache/version token ONCE here, while the source is
        # guaranteed present (the isfile guard above), and thread it into both
        # the thumb cache key and the X-Image-Version response header. Computing
        # it once avoids a second mtime syscall and keeps the disk-cache key
        # identical to the header the client already cached against.
        vtoken = image_version_token(ip)

        raw_w = 0; raw_h = 0
        bytes = if is_thumb
            raw_h, raw_w = raw_image_size(path)
            ensure_thumb_cached(db, id, path; token = vtoken)
        else
            full = ensure_full_cached(db, id, path; token = vtoken)
            raw_h, raw_w = full.raw_h, full.raw_w
            full.bytes
        end

        hdrs = ["Content-Type"    => "image/png",
                "Cache-Control"   => "private, max-age=31536000, immutable",
                "X-Image-Version" => vtoken,
                "X-Image-Width"   => string(raw_w),
                "X-Image-Height"  => string(raw_h)]
        HTTP.Response(200, hdrs, bytes)
    end

    @patch "/api/exposures/{id}/status" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)

        raw_status = get(body, :status, nothing)
        status = raw_status === nothing ? nothing : String(raw_status)

        if status !== nothing && status ∉ ("accepted", "rejected")
            return HTTP.Response(422,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "status must be 'accepted', 'rejected', or null")))
        end

        rows = sqlrows(DBInterface.execute(db,
            "SELECT id FROM exposures WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "image not found")))

        return with_idempotency(db, req) do
            DBInterface.execute(db,
                "UPDATE exposures SET status = ? WHERE id = ?", [status, id])

            result = apply_event!(InTransaction(), db, req;
                kind = "set_exposure_status",
                entity_type = "exposure", entity_id = id,
                payload = Dict(:status => status))
            _enqueue_broadcast_from_result!(result, "set_exposure_status", "exposure", id)

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:id => id, :status => status)))
        end
    end

    @post "/api/exposures/{id}/analyze" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = sqlrows(DBInterface.execute(db,
            "SELECT e.id, st.analysis_dir
             FROM exposures e JOIN samples s ON s.id = e.sample_id
             JOIN setups st ON st.id = e.setup_id
             WHERE e.id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "image not found")))
        analysis_dir = String(rows[1].analysis_dir)

        # Issue #9: route is now queue-shaped — `with_idempotency` honors the
        # `X-Client-Op-Id` header so a network retry returns the cached
        # response rather than running analyze a second time.
        # `analyze_exposure!` emits its `analyze_run` event with
        # `defer_broadcast=true`; we re-fetch that row below and enqueue a
        # post-commit broadcast carrying the spec'd `post_state` envelope so
        # foreign tabs converge without a refetch round-trip (issue #12).
        return with_idempotency(db, req) do
            # Pass `req` so analyze_exposure!'s durable analyze_run row carries
            # this request's client_op_id. We then filter by client_op_id
            # below — unambiguous, vs. the prior MAX(id)-sentinel approach
            # which could alias if a concurrent analyze interleaved between
            # the snapshot and the SELECT (review issue #15).
            # F-WIPE W1: the wipe/re-attach runs on this manual path too, so
            # capture the dropped-members report for the SSE post_state below.
            ares = analyze_exposure!(db, id, analysis_dir;
                                     defer_broadcast = true, req = req)

            new_hash = read_inputs_hash(db, id)
            client_op_id = get_client_op_id(req)
            # client_op_id is guaranteed non-nothing here because the queue
            # mutator always sends X-Client-Op-Id; a missing header would
            # have skipped the `with_idempotency` cache and we'd have no
            # disambiguation handle. Defend anyway — fall back to LIMIT 1
            # ORDER BY id DESC, which matches the prior behaviour.
            evt = if client_op_id === nothing
                sqlrows(DBInterface.execute(db,
                    """SELECT id, user_id, payload
                       FROM user_actions
                       WHERE action = 'analyze_run' AND entity_type = 'exposure'
                         AND entity_id = ?
                       ORDER BY id DESC LIMIT 1""", [id]))
            else
                sqlrows(DBInterface.execute(db,
                    """SELECT id, user_id, payload
                       FROM user_actions
                       WHERE action = 'analyze_run' AND entity_type = 'exposure'
                         AND entity_id = ? AND client_op_id = ?
                       LIMIT 1""", [id, client_op_id]))
            end
            if !isempty(evt)
                row = evt[1]
                payload_obj = ismissing(row.payload) ? nothing :
                              JSON3.read(String(row.payload)::String)
                # M0.4 mirror (see `_maybe_broadcast_event!`): a no-op reanalyze
                # (both skip flags true) must emit NO frame from this manual
                # route either — foreign tabs have nothing to converge on. This
                # hand-enqueued path previously bypassed the suppression that
                # the non-deferred `analyze_exposure!` path applies.
                noop = payload_obj !== nothing &&
                       get(payload_obj, :findpeaks_skipped, false) === true &&
                       get(payload_obj, :indexpeaks_skipped, false) === true
                if !noop
                    # F-WIPE W1: same post_state contract as the curation frames
                    # (shared builder, see `_enrich_curation_post_state`) — the
                    # assignment envelope carries the re-attached member ids, and
                    # `assignment_dropped` rides along only when the reanalyze
                    # dropped members. Without it, the frontend's envelope-absent-
                    # means-do-not-touch rule strands cached members on the
                    # cascaded (dead) index ids until a hard refetch.
                    post_state = _enrich_curation_post_state(db, id;
                        dropped_assignment_phases = ares.dropped_assignment_phases)
                    # The system-emitted row carries NULL client_id/op_id; pass
                    # the request's headers in the SSE frame so foreign tabs can
                    # still self-echo-filter (own client_id matches).
                    _enqueue_post_commit_broadcast!(
                        Int(row.id), "analyze_run", "exposure", id,
                        ismissing(row.user_id) ? nothing : Int(row.user_id),
                        get_client_id(req),
                        get_client_op_id(req),
                        ismissing(row.payload) ? nothing : String(row.payload);
                        post_state = post_state)
                end
            end

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(
                    :id                   => id,
                    :analyzed             => true,
                    :analysis_inputs_hash => new_hash,
                )))
        end
    end
end

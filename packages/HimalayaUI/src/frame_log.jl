const FRAME_LOG_DIR = ("analysis", "automatic_analysis", "sastool")
const FRAME_LOG_SETTLE_SECONDS = 60

const _FRAME_LOG_SAMPLE_SUMMARY = r"^Sample frames:\s*(\d+)\s+Rejected:\s*(\d+)\s*$"m
const _FRAME_LOG_ANY_SUMMARY = r"^(Sample|Buffer) frames:"m
const _FRAME_LOG_BUFFER_HEADER = "Variance of frames relative to the first frame for buffer files:"
const _FRAME_LOG_SAMPLE_HEADER = "Variance of frames relative to the first frame for sample files:"
const _FRAME_LOG_FRAME_LINE = r"^frame #\d+:"
const _FRAME_LOG_REJECTED_SUFFIX = "Frame rejected!!"

_frame_counts(total::Int, rejected::Int, final::Bool) =
    (total > 0 && 0 <= rejected <= total) ? (; total, rejected, final) : nothing

function _first_buffer_frame_counts(text::AbstractString)
    lines = split(text, '\n')
    start = findfirst(l -> startswith(strip(l), _FRAME_LOG_BUFFER_HEADER), lines)
    start === nothing && return nothing
    frames = String[]
    for l in @view lines[start+1:end]
        s = strip(l)
        occursin(_FRAME_LOG_FRAME_LINE, s) || break
        push!(frames, s)
    end
    _frame_counts(length(frames), count(endswith(_FRAME_LOG_REJECTED_SUFFIX), frames), false)
end

# physics: the beamline's averaged trace (`_tot.dat`) is the mean of the frames its software
# kept, so `total - rejected` is the number of frames the trace represents.
function parse_frame_log(text::AbstractString)
    m = match(_FRAME_LOG_SAMPLE_SUMMARY, text)
    m === nothing || return _frame_counts(parse(Int, m.captures[1]), parse(Int, m.captures[2]), true)
    has_summary = occursin(_FRAME_LOG_ANY_SUMMARY, text) || occursin(_FRAME_LOG_SAMPLE_HEADER, text)
    has_summary ? nothing : _first_buffer_frame_counts(text)
end

function frame_log_path(root_dir::AbstractString, image_pattern::AbstractString, stem::AbstractString)
    first_frame = first(splitext(replace(basename(String(image_pattern)), "{name}" => String(stem))))
    joinpath(String(root_dir), FRAME_LOG_DIR..., first_frame * ".log")
end

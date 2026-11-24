#!/usr/bin/env bash
set -euo pipefail

dist=${1:?dist path}
release=${2:?release path}
next_build=${3:?staged sysimage dir}
build=${4:?live sysimage dir}
service=${5:?systemd unit}
keep=${HIMALAYA_KEEP_RELEASES:-3}

releases=$(dirname "$release")
link_target="$(basename "$releases")/$(basename "$release")"

if [ ! -f "$release/index.html" ]; then
	echo "No frontend build at $release; nothing swapped." >&2
	exit 1
fi
if [ "$(cd "$(dirname "$releases")" && pwd -P)" != "$(cd "$(dirname "$dist")" && pwd -P)" ]; then
	echo "Releases must sit next to $dist; nothing swapped." >&2
	exit 1
fi

flip_dist() {
	ln -sfn "$1" "$dist.link"
	mv -T "$dist.link" "$dist"
}

if [ -d "$dist" ] && [ ! -L "$dist" ]; then
	legacy="00000000-000000-before-$(basename "$release")"
	if [ -e "$releases/$legacy" ]; then
		echo "$releases/$legacy already exists; nothing swapped." >&2
		exit 1
	fi
	mv -T "$dist" "$releases/$legacy"
	flip_dist "$(basename "$releases")/$legacy"
fi

previous=""
if [ -L "$dist" ]; then
	previous=$(readlink "$dist")
fi

flip_dist "$link_target"

if ! mv "$next_build/himalaya.so" "$build/himalaya.so"; then
	if [ -n "$previous" ]; then flip_dist "$previous"; else rm -f "$dist"; fi
	rm -rf "${release:?}"
	echo "Could not move the new sysimage into $build; the old frontend is back." >&2
	exit 1
fi
if ! mv "$next_build/julia_version" "$next_build/depot_path" "$build/"; then
	echo "The new sysimage is live, but its stamps did not move from $next_build to $build; make check-sysimage will report a mismatch." >&2
	exit 1
fi

live=$(basename "$link_target")
was_live=$(basename "${previous:-none}")
for dir in "$releases"/*/; do
	basename "$dir"
done | sort -r | tail -n +"$((keep + 1))" |
	while read -r old; do
		if [ "$old" != "$live" ] && [ "$old" != "$was_live" ]; then rm -rf "${releases:?}/$old"; fi
	done

echo
echo "Deployed. The running server keeps serving the old version until you restart it:"
echo
echo "    sudo systemctl restart $service"
echo

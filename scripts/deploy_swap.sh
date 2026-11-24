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
	mv "$dist" "$releases/00000000-000000"
	flip_dist "$(basename "$releases")/00000000-000000"
fi

previous=""
if [ -L "$dist" ]; then
	previous=$(readlink "$dist")
fi

flip_dist "$link_target"

if ! mv "$next_build/himalaya.so" "$build/himalaya.so"; then
	if [ -n "$previous" ]; then flip_dist "$previous"; else rm -f "$dist"; fi
	echo "Could not move the new sysimage into $build; the old frontend is back." >&2
	exit 1
fi
mv "$next_build/julia_version" "$next_build/depot_path" "$build/"

live=$(basename "$link_target")
find "$releases" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort -r | tail -n +"$((keep + 1))" |
	while read -r old; do
		[ "$old" = "$live" ] || rm -rf "${releases:?}/$old"
	done

echo
echo "Deployed. The running server keeps serving the old version until you restart it:"
echo
echo "    sudo systemctl restart $service"
echo

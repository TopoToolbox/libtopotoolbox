#!/bin/sh

git clone -c "core.sshCommand=ssh -i ~/.ssh/id_metrics -F /dev/null" --depth 1 git@github.com:TopoToolbox/metrics.git

echo "\"$(git rev-parse HEAD)\"" > timings
time -a -o timings -f "%e" cmake -B build/release -DCMAKE_BUILD_TYPE=Release -DTT_BUILD_TESTS=1 -DTT_DOWNLOAD_SNAPSHOTS=1
time -a -o timings -f "%e" cmake --build build/release
time -a -o timings -f "%e" ctest --test-dir build/release

jq -c -s '{project: "libtopotoolbox", timestamp: now, commit: .[0], metrics: {configure: .[1], build: .[2], test: .[3]}}' timings >> metrics/metrics.json

rm timings

cd metrics
git add metrics.json
git -c user.name=TT3MetricsBot -c user.email=metrics@example.com commit -m "Add libtopotoolbox metrics"
git push origin main

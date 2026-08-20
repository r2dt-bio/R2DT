# Local Traveler patches

Patches in this directory are applied on top of the pinned
`traveler_COMMIT` in [`../Dockerfile`](../Dockerfile), in the `traveler`
build stage, immediately after `git checkout`.

They exist because [cusbg/traveler](https://github.com/cusbg/traveler) is
currently dormant — the commit R2DT pins is also the head of upstream
`master`, which has had no new commits since November 2024. Patches here
are a stopgap, not a fork: each one corresponds to an open upstream PR and
should be deleted as soon as that PR is merged and `traveler_COMMIT` is
bumped past it.

## Current patches

| Patch | Upstream PR | Applies to | Remove when |
| --- | --- | --- | --- |
| `traveler-pr20-pk-per-bp.patch` | [cusbg/traveler#20](https://github.com/cusbg/traveler/pull/20) | `78c1cd63` | PR #20 is merged and `traveler_COMMIT` is bumped past it |

## Adding a patch

Generate it against the pinned commit so it applies from the Traveler
repository root:

```bash
git clone https://github.com/cusbg/traveler && cd traveler
git fetch origin pull/<N>/head:pr<N>
git diff <traveler_COMMIT>...pr<N> > traveler-pr<N>-<slug>.patch
```

Prepend a provenance header (upstream PR URL, commit, author, the commit
it applies to, and what it changes) — `git apply` ignores any text before
the first `diff --git` line. Then add a `COPY`-free one-liner to the
`traveler` stage's `git apply` invocation and a row to the table above.

## Verifying

Patches must apply cleanly to the pinned commit and compile in the same
image used by the build stage (`gcc:13`):

```bash
git -C traveler apply --check base_image/patches/<patch>
docker run --rm -v "$PWD/traveler:/t" -w /t/src gcc:13 make build
```

Traveler does not build on macOS/clang — `bits/stdc++.h` is a GCC-only
header — so verify inside the container, not on the host.

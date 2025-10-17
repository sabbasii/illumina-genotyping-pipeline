# Illumina Genotyping Pipeline

Linux-native, POSIX-compatible pipeline intended to run on WSL, native Linux, and macOS.

## Quick start
- Work inside a Linux filesystem path (e.g., `/home/<user>`). Avoid `/mnt/c/...` on WSL for performance.
- Raw data and outputs are **not** tracked in git (see `.gitignore`).
- Scripts and small config files *are* tracked.

## Environments
- WSL2 (Ubuntu), native Linux, macOS (zsh/bash)
- Dependencies will be listed in `docs/requirements.md` (to be added)

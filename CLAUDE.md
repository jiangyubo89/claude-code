# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Minimal Python/PyTorch utility scripts used alongside **idtracker.ai** (an animal tracking system). The environment runs Python 3.12.7 on Windows 11 with GPU support via CUDA.

## Running Scripts

```bash
python main.py   # Check CUDA availability
python jyb.py    # Same
```

## Environment

- **Python:** 3.12.7
- **Key dependencies (from log):** PyTorch, NumPy 2.0.1, OpenCV 4.10.0.84 (headless), PyQt6 6.10.2
- **idtracker.ai:** 6.0.13 — settings stored in Windows registry at `HKEY_CURRENT_USER\Software\idtrackerai`
- **Multiprocessing:** uses `spawn` method (Windows default)

## Notes

- No `requirements.txt` or virtual environment config present — dependencies are installed system-wide.
- No test framework configured.

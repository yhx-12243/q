#!/usr/bin/env python3

import hashlib
import requests
from argparse import ArgumentParser
from enum import Enum
from pathlib import Path
from shutil import rmtree
from subprocess import run

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--cargo-path', help='Cargo path of your system (default: ~/.cargo)', type=Path, default=Path.home() / '.cargo')
    parser.add_argument('--rustup-path', help='Rustup path of your system (default: ~/.rustup)', type=Path, default=Path.home() / '.rustup')
    parser.add_argument('--std-patch-server', help='Server to download std patch which updates frequently (default: https://43.138.56.99/rust-std-patches/)', default='https://43.138.56.99/rust-std-patches/')
    return parser.parse_args()

STD = ['core', 'alloc', 'std']

class PatchStatus(Enum):
    CLEAN = 0,
    PATCHED = 1,
    BROKEN = 2,
    DISCARD = 3,

def patch_inner(patch, path, is_std=False):
    assert patch.is_file()
    assert path.is_dir()

    hunks = []
    with open(patch) as f:
        for line in f.readlines():
            if line.startswith('index '):
                sH, tH = line.split()[1].split('..')
            elif line.startswith('--- '):
                sF = line[6:].rstrip()
            elif line.startswith('+++ '):
                tF = line[6:].rstrip()
                hunks.append((sF, tF, sH, tH))
    state = None

    def check_hunk(fn, ha):
        try:
            with open(path / fn) as f:
                ct = f.read()
        except FileNotFoundError:
            return ha.count('0') == len(ha)
        h = hashlib.sha1()
        ct = ct.encode()
        h.update(f'blob {len(ct)}'.encode())
        h.update(b'\0')
        h.update(ct)
        return h.hexdigest().startswith(ha)

    for sF, tF, sH, tH in hunks:
        match state:
            case None:
                match check_hunk(sF, sH), check_hunk(tF, tH):
                    case True, False:
                        state = PatchStatus.CLEAN
                    case False, True:
                        state = PatchStatus.PATCHED
                    case False, False:
                        state = PatchStatus.BROKEN
                        break
            case PatchStatus.CLEAN:
                if not check_hunk(sF, sH):
                    state = PatchStatus.BROKEN
                    break
            case PatchStatus.PATCHED:
                if not check_hunk(tF, tH):
                    state = PatchStatus.BROKEN
                    break
            case PatchStatus.BROKEN:
                break

    # remove std is not supported
    if state == PatchStatus.BROKEN:
        if is_std:
            state = PatchStatus.DISCARD
        else:
            rmtree(path)

    return patch, path, PatchStatus.PATCHED if state is None else state

def patch_std(identifier, patch, stdlib):
    return (
        f'\x1b[33m======== Applying \x1b[1;35m{identifier}\x1b[33m ========\x1b[0m',
        patch_inner(patch, stdlib / identifier, True),
    )

def patch_cargo(identifier, patch, crates_io):
    return (
        f'\x1b[36m======== Applying \x1b[1;35m{identifier}\x1b[36m ========\x1b[0m',
        patch_inner(patch, crates_io / identifier),
    )

def patch_git(identifier, patch, cargo_git):
    name, version = identifier.rsplit('-', 1)
    for d in cargo_git:
        if d.name.startswith(name):
            return (
                f'\x1b[32m======== Applying \x1b[1;35m{identifier}\x1b[32m ========\x1b[0m',
                patch_inner(patch, d / version),
            )
    raise FileNotFoundError(f'Cannot find {identifier}')

def main():
    args = parse_args()
    assert args.cargo_path.is_dir()
    assert args.rustup_path.is_dir()

    workspace = Path(__file__).parent
    assert workspace.is_dir()

    patches = workspace / 'patches'
    assert patches.is_dir()

    toolchain = run(['rustup', 'default'], capture_output=True, cwd=workspace)  \
        .stdout                                                                 \
        .split(None, 1)[0]                                                      \
        .decode()
    print(f'Use toolchain: \x1b[1;36m{toolchain}\x1b[0m\n')

    stdlib = args.rustup_path / 'toolchains' / toolchain / 'lib/rustlib/src/rust/library'
    assert stdlib.is_dir()

    crates_io = args.cargo_path / 'registry/src/index.crates.io-6f17d22bba15001f'
    assert crates_io.is_dir()

    cargo_git_dir = args.cargo_path / 'git/checkouts'
    cargo_git = []
    try:
        cargo_git = [d for d in cargo_git_dir.iterdir() if d.is_dir()]
    except:
        pass

    if args.std_patch_server:
        for std in STD:
            print(f'\x1b[35m======== Downloading \x1b[1;34m{std}\x1b[35m ========\x1b[0m\n')
            url = args.std_patch_server + ('' if args.std_patch_server.endswith('/') else '/') + std + '.patch'
            res = requests.get(url)
            with open(patches / f'{std}.patch', 'wb') as f:
                for chunk in res.iter_content(chunk_size=65536):
                    if chunk:
                        f.write(chunk)

    responses = []
    need_fetch = False
    for patch in patches.iterdir():
        name = patch.name
        if not name.endswith('.patch'):
            continue
        name = name[:-6]
        if name in STD:
            response = patch_std(name, patch, stdlib)
        elif (crates_io / name).is_dir():
            response = patch_cargo(name, patch, crates_io)
        else:
            response = patch_git(name, patch, cargo_git)
        if response:
            responses.append(response)
            if response[1][2] == PatchStatus.BROKEN:
                need_fetch = True

    if need_fetch:
        run(['cargo', 'fetch'], cwd=workspace)

    for (prompt, (patch, path, status)) in responses:
        print(prompt)
        match status:
            case PatchStatus.PATCHED:
                print('\x1b[1;32m√\x1b[0m')
            case PatchStatus.DISCARD:
                print(f'std package {path.name} is broken, discarded \x1b[1;33m(please use the NEWEST NIGHTLY version of Rust)\x1b[0m')
            case _:
                run(['git', '-C', str(path), 'apply', '--reject', str(patch)])
        print()

if __name__ == '__main__':
    main()

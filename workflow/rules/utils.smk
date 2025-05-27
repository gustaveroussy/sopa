def get_input_resolve(name: str, method_name: str):
    def _(wilcards):
        with getattr(checkpoints, f"patchify_{name}").get(**wilcards).output.patches_file.open() as f:
            resolve_paths = paths.temporary_boundaries_paths(f.read(), method_name)
            return [str(path.as_posix()) for path in resolve_paths]  # snakemake uses posix paths (fix issue #64)
    return _

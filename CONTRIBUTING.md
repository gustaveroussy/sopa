# Contributing to *sopa*

Contributions are welcome as we aim to continue improving the pipeline. For instance, you can contribute by:

- Opening an issue
- Discussing the current state of the code
- Making a Pull Request (PR)

If you want to open a PR, follow the following instructions.

## Installing `sopa` in editable mode

When contributing, installing `sopa` in editable mode is recommended. Also, we recommend installing the 'dev' extra.

For this, use one of the two following lines:

```sh
# with pip
pip install -e '.[dev]'

# or with uv
uv sync --all-extras --dev
```

## Coding guidelines

- Some code quality controls can be executed via pre-commit. You can install it via `pre-commit install`, and it will run the checks before any commit. You can also run `pre-commit` manually via `pre-commit run --all-files`
- Follow the [PEP8](https://peps.python.org/pep-0008/) style guide.
- Provide meaningful names to all your variables and functions.
- Document your functions and type your function inputs/outputs.
- Try as much as possible to follow the same coding style as the rest of the repository.

## Pull Requests

To add some new code to **sopa**, you should:

1. Fork the repository
2. Install `sopa` in editable mode with the 'dev' extra (see above)
3. Create your personal branch from `dev`
4. Implement your changes
5. Run tests via `pytest` (for coverage, use `pytest --cov --cov-config=pyproject.toml --cov-report=html` and open it with `open htmlcov/index.html`)
6. Run `pre-commit run --all-files` to ensure minimal code quality.
7. Commit and push changes
8. Create a pull request on the `dev` branch. Add explanations about your developed features, and wait for discussion and validation of your pull request

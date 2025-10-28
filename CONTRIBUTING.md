# Contributing to *sopa*

Contributions are welcome as we aim to continue improving the pipeline. For instance, you can contribute by:

- Opening an issue
- Discussing the current state of the code
- Making a Pull Request (PR)

If you want to open a PR, follow the following instructions.

## Installing `sopa` in editable mode

When contributing, installing `sopa` in editable mode is recommended. You'll need the `dev` extra to be able to run tests and pre-commits.

For this, use one of the two following lines (using `uv` is recommended, see [here](https://docs.astral.sh/uv/getting-started/installation/) how to install it):

```sh
# pip setup
pip install -e '.[dev]'

# uv setup
uv sync --all-extras --dev
```

## Coding guidelines

- Some code quality controls can be executed via pre-commit. You can install it via `pre-commit install`, and it will run the checks before any commit. You can also run it manually via `pre-commit run --all-files` to check the code quality at any time
- Follow the [PEP8](https://peps.python.org/pep-0008/) style guide.
- Provide meaningful names to all your variables and functions.
- Document your functions and type your function inputs/outputs.
- Try as much as possible to follow the same coding style as the rest of the repository.

## Pull Requests

To add some new code to **sopa**, you should:

1. Fork the repository
2. Install `sopa` in editable mode with the 'dev' extra (see above)
3. Create your personal branch from `main`
4. Implement your changes
5. Run tests via `uv run poe test_short`. Afterwards, you can open the test coverage report with `open htmlcov/index.html`. For the full tests, use `uv run poe test` (but it may take minutes).
6. Run `pre-commit run --all-files` to ensure minimal code quality.
7. Commit and push changes
8. Create a pull request on the `main` branch. Add explanations about your developed features, and wait for discussion and validation of your pull request

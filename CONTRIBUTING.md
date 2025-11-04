# Contributing to `sopa`

Contributions are welcome as we aim to continue improving the pipeline. For instance, you can contribute by:

- Opening an issue
- Discussing the current state of the code
- Making a Pull Request (PR)

If you want to open a PR, follow the following instructions.

## Installing `sopa` in editable mode

When contributing, installing `sopa` in editable mode is recommended. You'll need the `dev` extra to be able to run tests and pre-commits.

To do that, using `uv` is recommended. See [here](https://docs.astral.sh/uv/getting-started/installation/) how to install `uv`, and then install all the dependencies as below.

```sh
uv sync --all-extras --dev
```

## Coding guidelines

- Some code quality controls can be executed via pre-commit. You can install it via `pre-commit install`, and it will run the checks before any commit. You can also run it manually via `pre-commit run --all-files` to check the code quality at any time
- Follow the [PEP8](https://peps.python.org/pep-0008/) style guide.
- Provide meaningful names to all your variables and functions.
- Document your functions and type your function inputs/outputs.
- Try as much as possible to follow the same coding style as the rest of the repository.

## Pull Requests guidelines

To add some new code to `sopa`, you should:

### Branch setup
1. Fork the repository
2. Install `sopa` in editable mode with the 'dev' extra (see above)
3. Create your personal branch from `main`

### Code implementation and quality testing
4. Implement your changes
5. If you need to add documentation, add it under the `docs` directory. You can deploy it locally via `uv run poe docs`
6. Run tests via `uv run poe test_short`. You can open the test coverage report with `open htmlcov/index.html`. For the full tests, use `uv run poe test` (but it may take minutes).
7. If you fixed an issue or made an important change, mention it in the `CHANGELOG.md` file with your GitHub username.
8. Run `pre-commit run --all-files` to ensure minimal code quality.

### Create the PR
9.  Commit and push changes to your personal branch
10. Create a pull request on the `main` branch of the `sopa` repository. Add explanations about your developed features, and wait for discussion and validation of your pull request.

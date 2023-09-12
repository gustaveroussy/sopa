# Contributing to *sopa*

Contributions are welcome as we aim to continue improving the pipeline. For instance, you can contribute by:

- Reporting a bug
- Reporting some difficulties in using the pipeline
- Discussing the current state of the code
- Submitting a fix
- Proposing new features

## Contributing to the code, or debugging

Even though *sopa* is meant for flamingo, you can set up this repository locally for development or debugging purposes.

1. Install the dependencies.
2. Create your personal branch from `master`.
3. Make sure you read the coding guidelines below.
4. Implement your changes.
5. Create a pull request with explanations about your developed features. Then, wait for discussion and validation of your pull request.

## Coding guidelines

- Use the `black` formatter and `isort`. Their usage should be automatic as they are in the `pyproject.toml` file. Depending on your IDE, you can choose to format your code on save.
- Follow the [PEP8](https://peps.python.org/pep-0008/) style guide. In particular, use snake_case notations (and PascalCase for classes).
- Provide meaningful names to all your variables and functions.
- Document your functions and type your function inputs/outputs.
- Create your functions in the intended file, or create one if needed. See the project layout.
- Try as much as possible to follow the same coding style as the rest of the repository.
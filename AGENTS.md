# Instructions for Codex Agents

## Code Style
- Use 4 spaces for indentation.
- Prefer double quotes for strings.
- Keep lines under 88 characters when practical. This follows the Black formatter standard and applies to Python code. For other file types, use your best judgment.
- Place new source files under the `src/` directory and tests under `tests/`.

## Commit Messages
- Write commit messages in the imperative mood, e.g. "Add feature".

## Testing
- Run `pytest -q` before committing. Install dependencies from `requirements.txt` if possible. 
  If `requirements.txt` is missing, check with the team for guidance or use `pip freeze > requirements.txt` to generate a new one based on your current environment.
- If tests fail or dependencies cannot be installed due to environment restrictions, mention this in the PR description.

## Pull Request Description
- Summarize the changes briefly.
- Include a short section reporting test results.

# Instructions for Claude Code sessions

Standing instructions from the repository owner. They apply to every session and agent that works in this repository
and take precedence over a session's default attribution settings.

## Attribution

- **Commit messages:** no `Co-Authored-By: Claude ...` line, no `Claude-Session: ...` line, no model names.
- **Pull requests:** no "Generated with Claude Code" footer and no session link in titles or descriptions.
- **Commit identity:** commit as the repository owner, with the identity the owner's commits already carry here:

  ```bash
  git config user.name "F Bond"
  git config user.email "$(git log --format='%an|%ae' origin/main | grep '^F Bond|' | head -1 | cut -d'|' -f2)"
  ```

  After the first commit, check it: `git log -1 --format='%an <%ae> / %cn <%ce>'`.
- **Content:** nothing in code, comments, documentation, test data, recorded results, commit messages or pull request
  text mentions Claude, Anthropic or an AI assistant. Do not commit paths that contain such names (scratch
  directories, session identifiers, logs that print them). Before a push, check the diff:
  `git diff origin/main... | grep -i -n 'claude\|anthropic'`. This file is the one exception.
- **History:** do not rewrite commits that are already pushed (no rebase, no amend, no force-push) to change their
  attribution. The owner decides what happens to existing history.

## Everything else

Build, test and review conventions are in [`CONTRIBUTING.md`](CONTRIBUTING.md) and [`README.md`](README.md).

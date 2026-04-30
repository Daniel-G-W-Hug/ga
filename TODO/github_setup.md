# GitHub Repository Setup

## Branch Model

```text
main          ← releases only, strictly protected
  ↑
develop       ← owner's personal working branch, direct push allowed
  ↑
feature/xyz   ← coworker branches, merged into develop via pull request
```

## Branch Protection Rules

Configure in: **GitHub → Settings → Branches → Add rule**

### `main` — release branch

| Setting                               | Value                  | Reason                                 |
| ------------------------------------- | ---------------------- | -------------------------------------- |
| Require a pull request before merging | ✅ on                  | Nothing goes to `main` without a PR    |
| Required approvals                    | 1                      | Owner reviews before every release     |
| Require status checks to pass         | ✅ if CI is configured | Prevents broken releases               |
| Restrict who can push                 | ✅ owner only          | Coworkers cannot push to `main` at all |
| Allow force pushes                    | ❌ off                 | Release history must be immutable      |
| Allow deletions                       | ❌ off                 |                                        |

### `develop` — owner's working branch

| Setting                               | Value               | Reason                                                      |
| ------------------------------------- | ------------------- | ----------------------------------------------------------- |
| Require a pull request before merging | ✅ on               | Coworkers must submit a PR, not push directly               |
| Required approvals                    | 0 or 1              | Owner's choice — review coworker PRs before merging         |
| Restrict who can push                 | ✅ owner only       | Only owner pushes directly                                  |
| Add owner to bypass list              | ✅                  | Eliminates "bypassed rule violation" warning on direct push |
| Allow force pushes                    | ❌ off (for others) |                                                             |
| Allow deletions                       | ❌ off              |                                                             |

**Note on bypass list vs. push restriction:** "Restrict who can push" blocks everyone except
the named list. Adding the owner to the bypass list additionally makes direct pushes
explicit and intentional, so GitHub stops logging them as rule violations.

## Coworker Workflow

### Starting new work

```bash
# Get the latest state of develop
git checkout develop
git pull origin develop

# Create a feature branch from develop
git checkout -b feature/my-feature
```

### Submitting work

```bash
# Push the feature branch to GitHub
git push origin feature/my-feature

# Then open a Pull Request on GitHub:
# feature/my-feature → develop
```

The owner reviews the PR and merges it into `develop`.

### Keeping a feature branch up to date

```bash
# While working on a feature branch, incorporate new changes from develop
git fetch origin
git rebase origin/develop
# or: git merge origin/develop
```

### Branch naming conventions (suggestions)

| Prefix      | Purpose            | Example                        |
| ----------- | ------------------ | ------------------------------ |
| `feature/`  | New functionality  | `feature/vtk-3d-visualization` |
| `fix/`      | Bug fix            | `fix/rotor-normalization`      |
| `refactor/` | Code restructuring | `refactor/split-active-common` |
| `test/`     | Test additions     | `test/pga3d-congruence`        |

## Owner Workflow

```bash
# Push directly to develop — no PR needed
git push origin develop

# Release: when develop is stable, open a PR on GitHub
# develop → main
# Review, then merge
```

## Summary

| Who      | Branch        | How                       |
| -------- | ------------- | ------------------------- |
| Owner    | `develop`     | Direct push               |
| Owner    | `main`        | PR from `develop`         |
| Coworker | `feature/xyz` | Direct push to own branch |
| Coworker | `develop`     | PR only — no direct push  |
| Coworker | `main`        | No access                 |

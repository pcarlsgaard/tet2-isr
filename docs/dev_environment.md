# Working with this project in Docker + VS Code

A checklist for the setup this project was built in: a VS Code dev container whose files live in a Docker
**named volume**. Read "How it is set up" once, then use the numbered procedures.

## How it is set up

| Thing | Where it lives | Survives |
|---|---|---|
| Project files (`/workspaces/typescript-node/tet2-isr`) | Docker volume `typescript-node` inside Docker Desktop's Linux VM. **Not visible in Finder.** | stop, restart, rebuild, remove container |
| Installed tools (`gh`, apt packages) | the container's own layer | stop and restart only. **Lost on rebuild or remove.** |
| `data/raw/`, `data/processed/`, `.venv/` | inside the project folder (volume) but gitignored | same as project files; recreate with `make public` / `make setup` |
| Code, docs, results tables, history | git | everything, once pushed to GitHub |

Rule: **GitHub is the source of truth.** The volume is a working copy. Anything not pushed can be lost if the volume is
deleted (Docker Desktop → Volumes → delete, or `docker volume rm typescript-node`).

## 1. First time: sign in to GitHub from the container

Run in the VS Code terminal inside the container:

```bash
gh --version                      # if "command not found", see "Reinstall gh" below
gh auth login                     # choose: GitHub.com → HTTPS → Login with a web browser
gh auth setup-git                 # makes plain `git push` use that login
git config --global user.name  "Your Name"
git config --global user.email "you@example.com"
```

Try `git push` before logging in; VS Code often forwards your Mac's Git credentials automatically.

Reinstall `gh` after a container rebuild:

```bash
sudo mkdir -p -m 755 /etc/apt/keyrings
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | sudo tee /etc/apt/keyrings/githubcli-archive-keyring.gpg >/dev/null
echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | sudo tee /etc/apt/sources.list.d/github-cli.list >/dev/null
sudo apt-get update && sudo apt-get install -y gh
```

## 2. Save your work to GitHub (do this often, and always before closing)

```bash
cd /workspaces/typescript-node/tet2-isr
git status                        # what changed
git add -A
git commit -m "describe what you did"
git push origin HEAD                          # pushes the branch you are on
```

Check it worked: open https://github.com/pcarlsgaard/tet2-isr and look for your commit.

`main` is the default branch and is protected (pull request required, no force-push, no deletion). Do your work on a
branch and merge through a pull request:

```bash
git checkout main && git pull origin main
git checkout -b my-change                     # one branch per task
# ...work, commit...
git push -u origin my-change
gh pr create --fill                           # then review and merge on GitHub, or: gh pr merge --squash
git checkout main && git pull origin main
```

The repo is public: do not commit unpublished or patient-derived data, tokens, or anything you would not post.

If a push is rejected, run `git pull --rebase origin main` first, then push again.

## 3. Close everything and free memory

Docker Desktop's background VM holds memory while it runs, even if the container is idle.

1. Save to GitHub (section 2). Confirm `git status` shows nothing uncommitted and `git log origin/main..HEAD` is empty (or your branch is pushed).
2. VS Code: **File → Close Remote Connection** (this stops the container).
3. Menu bar whale icon → **Quit Docker Desktop**. Memory is released at once. The volume is kept.

Optional memory cap: Docker Desktop → **Settings → Resources → Memory** → set 4–6 GB → Apply & Restart.
The full pipeline fits in that.

## 4. Reopen later

1. Start **Docker Desktop** and wait until it says running.
2. Open VS Code → **File → Open Recent** → pick the `typescript-node [Dev Container]` entry
   (or Command Palette → *Dev Containers: Reopen in Container*).
3. Terminal:
   ```bash
   cd /workspaces/typescript-node/tet2-isr
   git pull origin main          # in case you worked from the Mac
   make test                     # 23 tests, a few seconds
   ```
   If `make test` says the venv is missing: `make setup`.

## 5. Keep a copy on your Mac too (GitHub is the hub)

The container cannot write to a Mac folder directly, so sync both through GitHub.

On the **Mac** terminal (once):

```bash
git clone https://github.com/pcarlsgaard/tet2-isr.git ~/tet2-isr
cd ~/tet2-isr
```

Then, from either place:

```bash
git pull origin main              # before you start
# ...work...
git add -A && git commit -m "..." && git push origin HEAD
```

Pull before starting, and avoid editing the same branch in both places at the same time.
`data/` and `.venv/` are not in git; on the Mac run `make setup` (needs Python ≥ 3.10) and `make public` to recreate them.

## 6. Backups that do not depend on GitHub

**One-file git backup** (small, about 20 MB). In the container:

```bash
cd /workspaces/typescript-node/tet2-isr
git bundle create ../tet2-isr-backup.bundle --all
```

Copy it to the Mac (Mac terminal; find the container name with `docker ps`):

```bash
docker cp <container-name>:/workspaces/typescript-node/tet2-isr-backup.bundle ~/Desktop/
```

Restore later: `git clone ~/Desktop/tet2-isr-backup.bundle tet2-isr`.

**Whole-volume backup** (includes the 1.4 GB of downloaded data; Mac terminal):

```bash
docker run --rm -v typescript-node:/v -v "$PWD":/b alpine tar czf /b/typescript-node-backup.tgz -C /v .
```

**Single files in VS Code:** right-click in the Explorer → *Download*. Skip `data/` and `.venv/`.

## 7. Safety habits for agent work

- Let the agent work on a branch; review the diff before merging into `main`; turn on branch protection for `main` on GitHub.
- Use a fine-grained GitHub token limited to this repo with an expiry if you do not want a full login:
  `gh auth login --with-token < token.txt`, then delete the file. Revoke the token when you pause the project.
- Do not mount your Mac home folder, SSH keys or cloud credentials into the container.
- The container currently has open outbound network access and passwordless `sudo`. For stricter isolation, use a
  Python dev container with a firewall allowlist (GitHub, PyPI, NCBI, Ensembl, MGI, Broad, Anthropic API) and a
  non-root user. Not done yet.
- Put unpublished or patient-derived data in a read-only mount, not in an environment with open network access.

## 8. Common problems

| Symptom | Fix |
|---|---|
| `git push` asks for a username / "terminal prompts disabled" | do section 1 |
| `gh: command not found` | container was rebuilt; reinstall (section 1) |
| `make: .venv/bin/python: not found` | `make setup` |
| pipeline cannot find `data/raw/...` | `make public` (re-downloads about 1.4 GB, runs in about a minute after download) |
| files missing after a rebuild | check you opened the same volume: `docker volume ls` should list `typescript-node` |
| Docker still uses lots of memory | quit Docker Desktop; set the memory cap (section 3) |

## 9. Where to look

`README.md` (status and results), `reports/research_report.html` (real-data report), `reports/executive_summary.md`,
`docs/data_audit.md`, `docs/assumptions_and_limits.md`, `docs/decisions.md`, `AGENTS.md`.

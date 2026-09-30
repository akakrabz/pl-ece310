# ECE 329 notes site — how to build, preview, and host

This folder is a **Quartz 5** site (`quartz.config.yaml`, `quartz.ts`, `content/`, `quartz/static/demos/`, `quartz/styles/custom.scss`).
The notes are plain Markdown with Obsidian-style wikilinks and callouts, so the `content/` folder also opens directly as an Obsidian vault.

## 1. Build and preview (first time)

Requirements: Node ≥ 22, npm ≥ 10.9 (you have Node 22.23), internet access for the two install steps.

```bash
cd "~/Nextcloud/Notes/ECE 329/quartz"
npm ci                        # Quartz's own dependencies (includes the @quartz-community npm plugins)
npm run install-plugins       # generates .quartz/plugins/index.ts from quartz.config.yaml — required once, and after editing the plugin list
npx quartz build --serve      # http://localhost:8080  (rebuilds on save)
```

`npx quartz build` alone writes the static site to `public/`.

Note: `npx quartz plugin install` is only for plugins fetched from Git; it does *not* regenerate the index for npm plugins
(it just says "No quartz.lock.json found"). If the build fails with `Could not resolve "./.quartz/plugins"`, run `npm run install-plugins`.

## 2. Before hosting

`configuration.baseUrl` in `quartz.config.yaml` is set to `pl-ece329.vops.ch` (host plus any sub-path, no protocol, no trailing slash). It only affects the sitemap and social-preview URLs.

If the site lives under a sub-path (e.g. `example.com/ece329/`), build with `npx quartz build --baseDir /ece329`.

## 3. Hosting: the one rule your server needs

Quartz emits `page.html` files but links to them **without** the extension (`/concepts/gauss-law`), and folders are served with a trailing slash. Your web server must try `$uri`, then `$uri.html`, then `$uri/`:

```nginx
# nginx
location / {
    try_files $uri $uri.html $uri/ =404;
}
error_page 404 /404.html;
```

```caddyfile
# Caddy
try_files {path} {path}.html {path}/ =404
```

Apache: `Options -MultiViews` plus `RewriteEngine On; RewriteCond %{REQUEST_FILENAME}.html -f; RewriteRule ^(.*)$ $1.html [L]`.

Copy `public/` to the web root (or the sub-path). Opening `public/index.html` directly from disk (file://) does not work — search, graph and page previews fetch JSON over HTTP.

### Docker / Portainer (the current hosting)

`docker-compose.yml` in this folder is a self-contained stack for Docker standalone (Proxmox → Portainer): a `builder` container (`node:22`) clones `https://github.com/akakrabz/quartz-329`, checks it every 5 minutes and rebuilds the site into a shared volume whenever `master` moves; an nginx container serves that volume on port **8329** with the `try_files` rule above. Builds go to a staging directory and are renamed into place, so the old site keeps serving during a build and after a failed one.

- Deploy: Portainer › Stacks › Add stack › Web editor, paste the file (or `docker compose up -d`). Point the reverse proxy for `pl-ece329.vops.ch` at port 8329. No secrets — the repo is public.
- Publish a change: `git push`. The site follows within `SYNC_INTERVAL` (first start ≈ 2 min for `npm ci`; content-only rebuilds ≈ 15 s; `npm ci` only reruns when `package.json`/`package-lock.json` changed).
- Logs (builder): one `[sync …]` line per event — `new commit … building`, `site updated to …`, or `BUILD FAILED …` followed by the build error. A broken commit is not retried until the next push, so fix and push again.
- To point at a different repository, change `REPO_URL` and delete the `repo` volume; a different `BRANCH` needs no volume reset.

External requests the site makes at page load: Google Fonts (theme fonts) and jsdelivr (KaTeX CSS, loaded by the latex plugin). Analytics are disabled. To go fully self-hosted later, set `theme.fontOrigin: local` and vendor the KaTeX CSS.

## 4. Writing conventions (so new pages match)

- One page per lecture in `content/<unit>/NN-slug.md`; frontmatter `title`, `description`, `tags`, `lecture`.
- Concept pages in `content/concepts/`, problems in `content/problems/`, demos in `content/demos/`.
- Link with full paths: `[[concepts/gauss-law|Gauss's law]]`, `[[1-electrostatics/03-gauss-law-at-work#3-the-three-symmetries|Lecture 3 §3]]`. Inside tables escape the pipe: `[[page\|text]]`.
- Math: `$…$` inline, `$$…$$` on its own lines (every line prefixed with `> ` inside a callout). Use `\lvert x\rvert` instead of `|x|` inside tables. No custom macros — the vault must also render in Obsidian.
- Callouts: the standard Obsidian types plus this site's own `key`, `recipe`, `trap`, `exam`, `intuition`, `derivation` (styled in `quartz/styles/custom.scss`). Append `-` to the type to fold by default.
- Figures: inline `<figure class="ece-fig">…SVG…</figure>` blocks with **no blank lines inside** and **a blank line after** (an HTML block runs until a blank line — without one, the next callout or equation is swallowed as raw HTML and the build fails); strokes use `currentColor` and the CSS variables `--accent`, `--accent2`, `--hi`, `--muted` so they follow dark mode. The generator for the existing figures is in `tools/figs.py`.
- Demos: standalone HTML in `quartz/static/demos/<name>/index.html`, embedded with `<iframe src="/static/demos/<name>/">` inside `<div class="ece-demo">`.
- Explorer order comes from file names (numeric prefixes), see `quartz.ts`; titles stay clean.

## 5. Checking a page without building

`tools/check.py` (Python 3, needs PyYAML and the local KaTeX copy path set at the top) validates every wikilink and heading anchor and compiles every equation with KaTeX in strict mode — the same checks run before this delivery. `tools/assemble.py` inlines figures from `tools/figs/` into pages that contain `<!-- fig:name -->` markers, if you keep the figures separate.

## 6. What's here (build 2, 2026-09-29 — Unit 1 through the Exam 1 scope)

- Home page with the course map and conventions.
- Toolkit: coordinates & differential elements, vector-calculus cheat sheet, units & constants, errata in the course materials.
- Unit 1, Lectures 1–10 written in full (everything Exam 1 covers); Lecture 11 and Units 2–4 outlined on their index pages.
- 25 concept pages (the graph's hubs): the original 16 plus potential, boundary conditions, Poisson's equation, conductors, polarization, permittivity, capacitance, conductance, electrostatic energy.
- 5 worked problems: one per FA26 Exam 1 problem (#1–#4, all re-parameterized) plus the Lecture 3 flux challenge.
- 1 interactive demo (point charges + Gaussian loop).
- 25 original figures (`tools/figs/`, generated by `tools/figs.py`).

Every page passed `tools/check.py` (53 pages, 422 wikilinks/anchors, 3395 KaTeX expressions) and an independent physics review before delivery.

The course's own canvas demos (`Suppliment/Websites/*.html`, `smithchart.html`) were **not** copied into the public tree: they are saved from the course's login-only site and carry no license. They can be dropped into `quartz/static/demos/course-apps/` for a private build.

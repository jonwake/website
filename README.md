# Jon Wakefield — MkDocs Site

This repository is a [MkDocs](https://www.mkdocs.org/) version of Jon Wakefield's faculty website, built with the [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) theme.

## Site structure

```
.
├── mkdocs.yml            # MkDocs configuration (theme, nav, plugins)
├── requirements.txt      # Python dependencies
├── docs/                 # Markdown source for each page
│   ├── index.md
│   ├── research.md
│   ├── people.md
│   ├── software.md
│   ├── teaching.md
│   ├── spatial-demography.md
│   ├── regression-methods.md
│   ├── personal.md
│   ├── images/           # Photos used on pages
│   └── files/            # PDFs (lecture slides, COVID method papers, etc.)
└── .github/workflows/deploy.yml   # Auto-deploy to GitHub Pages on push
```

## Local development

Create a virtual environment, install dependencies, and serve the site locally with live reload:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
mkdocs serve
```

Then open <http://127.0.0.1:8000> in your browser. Edits to any `docs/*.md` file or to `mkdocs.yml` will reload automatically.

## Building a static site

```bash
mkdocs build
```

The output is written to `site/` (gitignored).

## Deployment

A GitHub Actions workflow (`.github/workflows/deploy.yml`) builds the site on every push to `main`/`master` and publishes it to GitHub Pages. After pushing to GitHub:

1. Open the repo's **Settings → Pages**.
2. Under "Build and deployment", set **Source** to **GitHub Actions**.
3. Push a commit (or run the `Deploy MkDocs to GitHub Pages` workflow manually); the site will be available at `https://<your-org-or-user>.github.io/<repo-name>/`.

Remember to update `site_url` in `mkdocs.yml` once you know the final URL.

## Adding or editing content

- Edit the `.md` file under `docs/` for the page you want to change.
- To add a new page, create a new `.md` file in `docs/` and add an entry to the `nav:` section of `mkdocs.yml`.
- Drop new images into `docs/images/` and reference them with relative paths like `![Alt text](images/myphoto.jpg)`.
- Drop new PDFs into `docs/files/` and link them with `[Slides](files/my-slides.pdf)`.

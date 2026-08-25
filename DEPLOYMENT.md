# 🌐 Deploy HEMOX with GitHub Pages

## 1. Verify the repository root

After extracting the upload package, upload its **contents**, not the ZIP itself and not an outer `HEMOX` folder.

The **Code** page of `Ahmed-Reda-SILO/HEMOX` must show this structure directly:

```text
.github/
.nojekyll
CHANGELOG.md
CITATION.cff
DEPLOYMENT.md
README.md
index.html
```

If GitHub instead shows `HEMOX/` or `HEMOX-v1.7.12-repository.zip`, the files are one level too deep and Pages will return 404.

## 2. Enable GitHub Actions as the Pages source

1. Open the `HEMOX` repository.
2. Select **Settings → Pages**.
3. Under **Build and deployment → Source**, select **GitHub Actions**.
4. Open **Actions → Deploy HEMOX to GitHub Pages**.
5. Select **Run workflow**, choose `main`, and run it.
6. Wait until the workflow and `github-pages` deployment both show a green check mark.

The published address is:

<https://ahmed-reda-silo.github.io/HEMOX/>

## 3. If the workflow does not appear

The workflow file is not at the required path. Confirm that it is exactly:

```text
.github/workflows/pages.yml
```

It must not be located at `HEMOX/.github/workflows/pages.yml` inside the repository.

## 4. Simple branch-based alternative

Because HEMOX is a static single-file application, it can also be published without a workflow:

1. Confirm that `index.html` is at the repository root on `main`.
2. Open **Settings → Pages**.
3. Select **Deploy from a branch**.
4. Choose `main` and `/(root)`.
5. Select **Save** and wait for the deployment to finish.

Use only one publishing method at a time. The GitHub Actions method is recommended for this repository.

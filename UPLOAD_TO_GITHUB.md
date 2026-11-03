# Upload this repository

The folder `atmosphere-contrast` is the repository root. Keep its internal
structure. Upload the extracted files, rather than the delivery ZIP.

## Option A: GitHub website

1. Create a repository named `atmosphere-contrast` in your GitHub account.
2. Do not generate another README, license, or gitignore; they are included.
3. Choose **uploading an existing file**, or **Add file > Upload files**.
4. Open the extracted `atmosphere-contrast` folder on your Mac. Drag its
   contents into the upload area so `README.md` is at the repository root.
5. In Finder, press **Command + Shift + .** to display dotfiles. Include
   `.github`, `.gitignore`, `.gitattributes`, and `.editorconfig` as well.
6. Use the commit message `Add documented atmosphere contrast review module`.
7. Commit the files, then open the repository's **Actions** tab for test
   results. The included workflow is configured to run after an upload; a
   workflow in this ZIP has not already run on GitHub.

If the browser does not upload dotfiles or the workflow folder, use Option B.
Gitignore rules are applied by Git; a browser upload still requires selecting
only the intended repository files. This supplied folder contains no model
downloads or private results.

## Option B: GitHub Desktop

1. Extract the ZIP and keep the `atmosphere-contrast` folder in a permanent
   location on your Mac.
2. In GitHub Desktop, choose **File > Add Local Repository** and select it.
3. If told that it is not a Git repository, use the offered **create a
   repository here** link. Keep this existing folder as the repository root.
4. Keep the existing README, MIT license, and gitignore.
5. Review the changes and commit them with the message above.
6. Click **Publish repository**, select your account and visibility, and
   publish. Share the resulting repository URL with Vivek.

## Included license and scope

An MIT license is included for the repository code/documentation, naming
Amanpreet Singh. External model data and transmission curves are not bundled
or relicensed. The README describes this as a review module, not the complete
Atmosphere Lab app or a completed upstream pyNOMIC integration.

Official upload instructions:
<https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository>


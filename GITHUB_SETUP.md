# Prepare GitHub hosting

1. Set the real maintainer name/email in DESCRIPTION and choose the LICENSE.
2. Install dependencies, run package tests/checks and review VALIDATION.md.
3. Run the supplied biological reference script with the external files.
4. Create an empty GitHub repository named LiDAR under the desired owner.
5. From the unpacked repository:

```sh
git remote add origin https://github.com/OWNER/LiDAR.git
git push -u origin main
```

The local Git repository is initialized with a preparation commit. Its commit
identity is a neutral packaging identity, not a claim of scientific authorship.
No GitHub repository has been created remotely and no files have been pushed.
Do not add the input QS or the >100 MB ligand-target RDS to ordinary Git.
External resource hosting and license/source metadata can be added after the
owner supplies the approved distribution location and database provenance.

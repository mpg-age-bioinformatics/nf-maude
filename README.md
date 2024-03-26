# nf-maude

Run the workflow:
```
nextflow run mpg-age-bioinformatics/nf-maude -r 1.0.0 -params-file params.json
```
or
```
git clone https://github.com/mpg-age-bioinformatics/nf-maude.git
nextflow run nf-maude -params-file params.json
```

## Contributing

Make a commit, check the last tag, add a new one, push it and make a release:
```
git add -A . && git commit -m "<message>" && git push
git describe --abbrev=0 --tags
git tag -e -a <tag> HEAD
git push origin --tags
gh release create <tag>
```

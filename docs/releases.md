# Release process

- All R2DT releases are [available](https://github.com/RNAcentral/R2DT/releases) on GitHub.
- R2DT uses [git flow](https://github.com/RNAcentral/R2DT/wiki) for managing the release process.
- See an example [GitHub issue](https://github.com/RNAcentral/R2DT/issues/86) listing the steps required for a release.

## git workflow

R2DT uses a popular git workflow that's often called
`git flow`. See the [2010 blog post by Vincent
Driessen](http://nvie.com/posts/a-successful-git-branching-model/)
that describes it. We will use it with the difference that we
don't mind having feature branches on `origin`.

In what follows, first we'll give concise-ish examples of the flow for normal
development, making a release, and making a "hotfix". A summary of the
principles and rationale follows the examples.

### Normal development

Generally, for any changes you make to our code, you will make on a
*feature* branch, off of `develop`. So first you create your branch:

```bash
   $ git checkout -b myfeature develop
```

Now you work, for however long it takes. You can make commits on your
`myfeature` branch locally, and/or you can push your branch up to the
origin and commit there too, as you see fit.

When you're done, and you've tested your new feature, you merge it to
`develop` (using `--no-ff`, which makes sure a clean new commit object
gets created), and delete your feature branch:

```bash
   $ git checkout develop
   $ git merge --no-ff -m "Merges myfeature branch into develop" myfeature
   $ git branch -d myfeature
   $ git push origin --delete myfeature
   $ git push origin develop
```

#### Small features: single commits can be made to `develop`

Alternatively, if you're sure your change is going to be a single
commit, you can work directly on the `develop` branch.

```bash
   $ git checkout develop
     # make your changes
   $ git commit
   $ git push origin develop
```

#### Big features: keeping up to date with `develop`

If your work on a feature is taking a long time (days, weeks...), and
if the `develop` trunk is accumulating changes you want, you might
want to periodically merge them in:

```bash
   $ git checkout myfeature
   $ git merge --no-ff -m "Merges develop branch into myfeature" develop
```

________________________________________________________________
### **Making a release**

To make a release, you're going to make a *release branch* of the
code, and of any other repos it depends on (currently none, but we
may update this in the future if we want to couple R2DT development
with active development of other packages it uses, e.g. Traveler or Ribovore).
You assign appropriate version numbers to the appropriate files in the
release branch, test and stabilize. When everything is ready, you merge to `master` and tag
that commit with the version number; then you also merge back to
`develop`, and delete the release branch.

For example, here's the git flow for an R2DT release.
Suppose R2DT is currently at v1.1.5 and we decide this release will
be R2DT 1.2. We first make a new release from R2DT's `develop` branch:

```bash
   $ cd r2dt
   $ git checkout develop # only necessary if you're not already on develop
   $ git checkout -b release-1.2 develop
     # change version number and dates in Readme.md and utils/shared.py
   $ git commit -a -m "Version number bumped to 1.2"
     # do and commit any other work needed to test/stabilize R2DT release.
```

Then merge the release branch as follows:

```bash
   $ git checkout master
   $ git merge --no-ff -m "Merges release-1.2 branch into master" release-1.2
   # Now merge release branch back to develop...
   $ git checkout develop
   $ git merge --no-ff -m "Merges release-1.2 branch into develop" release-1.2
   $ git push
   $ git branch -d release-1.2
   $ git push origin --delete release-1.2
```

Create a release using GitHub web interface that allows to include release notes.

________________________________________________________________
### **Fixing bugs: "hotfix" branches**

If you need to fix a critical bug and make a new release immediately,
you create a `hotfix` release with an updated version number, and the
hotfix release is named accordingly: for example, if we screwed up
R2DT 1.2, `hotfix-1.2.1` is the updated 1.2.1 release.

A hotfix branch comes off `master`, but otherwise is much like a
release branch.

```bash
   $ cd r2dt
   $ git checkout -b hotfix-1.2.1 master
     # bump version number to 1.2.1; also dates, copyrights
   $ git commit -a -m "Version number bumped to 1.2.1"
```

Now you fix the bug(s), in one or more commits. When you're done, the
finishing procedure is just like a release:

```bash
    $ git checkout master
    $ git merge --no-ff -m "Merges hotfix-1.2.1 branch into master" hotfix-1.2.1
    $ git tag -a r2dt-1.2.1
    $ git checkout develop
    $ git merge --no-ff -m "Merges hotfix-1.2.1 branch into develop" hotfix-1.2.1
    $ git push
    $ git branch -d hotfix-1.2.1
    $ git push origin --delete release-1.2
```


________________________________________________________________
### **Summary of main principles**

There are two long-lived R2DT branches: `origin/master`, and `origin/develop`. All other branches
have limited lifetimes.

`master` is stable.  Every commit object on `master` is a tagged
release, and vice versa.

`develop` is for ongoing development destined to be in the next
release. `develop` should be in a _close-to-release_ state.

We make a *feature branch* off `develop` for any nontrivial new work --
anything that you aren't sure will be a single commit on `develop`. A
feature branch:
 * comes from `develop`
 * is named anything informative (except `master`, `develop`, `hotfix-*` or `release-*`)
 * is merged back to `develop` (and deleted) when you're done
 * is deleted once merged

We make a *release branch* off `develop` when we're making a release.
A release branch:
 * comes from `develop`
 * is named `release-<version>`, such as `release-1.2`
 * first commit on the hotfix branch consists of bumping version/date/copyright
 * is merged to `master` when you're done, and that new commit gets
   tagged as a release
 * is then merged back to `develop` too
 * is deleted once merged

We make a *hotfix branch* off `master` for a critical immediate fix to
the current release. A hotfix branch:
 * comes from `master`
 * is named `hotfix-<version>`, such as `hotfix-1.2.1`
 * first commit on the hotfix branch consists of bumping version/date/copyright
 * is merged back to `master` when you're done, and that new commit
   object gets tagged as a release.
 * is then merged back to `develop` too
 * is deleted once merged

#### Much of the text above was borrowed and modified with permission from Sean Eddy, from the [HMMER GitHub repository](https://github.com/EddyRivasLab/hmmer/wiki/Git-workflow).

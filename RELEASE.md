# Release Instructions

Or, more directly, how to update the PyPI file.


## 1. Bump the version number in pyproject.toml

This is important and the PyPI upload will fail if you don't do this. You can update the this file either on GitHub or you can push the changes from your local machine. 

This package uses semantic versioning. From py-pkgs.org:

- Patch release (0.1.0 -> 0.1.1): patch releases are typically used for bug fixes, which are backward compatible. Backward compatibility refers to the compatibility of your package with previous versions of itself. For example, if a user was using v0.1.0 of your package, they should be able to upgrade to v0.1.1 and have any code they previously wrote still work. It’s fine to have so many patch releases that you need to use two digits (e.g., 0.1.27).

- Minor release (0.1.0 -> 0.2.0): a minor release typically includes larger bug fixes or new features that are backward compatible, for example, the addition of a new function. It’s fine to have so many minor releases that you need to use two digits (e.g., 0.13.0).

- Major release (0.1.0 -> 1.0.0): release 1.0.0 is typically used for the first stable release of your package. After that, major releases are made for changes that are not backward compatible and may affect many users. Changes that are not backward compatible are called “breaking changes”. For example, changing the name of one of the modules in your package would be a breaking change; if users upgraded to your new package, any code they’d written using the old module name would no longer work, and they would have to change it.

## 2. Create a release on GitHub

Navigate to the repository on GitHub and fine the 'Releases' section in the middle left of the repostories home page. 


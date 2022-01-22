# 1 Installation
***

## 1.1 Building the stable package
Note that at the early stages, stable package isn't bug free - tracking issues and fixing bugs will still happen. The most recent stable build is 0.0.2-stable

To install, simply run
```
$ pip3 install stable_builds/flares-segmentation-tlincke125-0.0.2-stable.tar.gz
```

## 1.2 Building the latest package 
First, checkout the **master** branch. Then, in the root directory, build the package using the following command:

```
$ python3 -m build
``` 

This should have produced two files in the ./dist folder.

To install the package, go into the dist folder and install the package using pip:

```
$ cd dist
$ pip3 install flares-segmentation-tlincke125-<version>.tar.gz
```

Where \<version\> is the highest version (or any other version you wish to install).

## 1.3 Building documentation
Documentation is generated using pdoc. To build the documentation, run:

```
$ pdoc -f --html -o ./docs/documentation ./src/flares -c latex_math=True
```

Then open docs/documentation/flares/index.html in your prefered browser

# 2 Getting Started
***
See ./docs/guides/getting_started.ipynb

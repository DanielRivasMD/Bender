<!-- TODO: update documation -->

## Bender is a set of Genomic tools for deployment on Slurm managed systems

![](assets/Bender.png)
![](assets/Slurm.png)


## Table of contents

  - [Table of contents](#table-of-contents)
  - [Overview](#overview)
  - [Installation](#installation)
    - [Via Homebrew (for macOS)](#via-homebrew-for-macos)
    - [Via APT (for Debian-based Linux distros)](#via-apt-for-debian-based-linux-distros)
    - [From Github release](#from-github-release)
  - [Documentation](#documentation)
    - [Usage](#usage)
  - [Examples](#examples)
    - [Example `bender` config](#example-bender-config)
  - [_Bender_ for the curious](#bender-for-the-curious)
  - [Acknowledgements](#acknowledgements)
  - [License](#license)


## Overview


## Installation


<!-- TODO: -->
**Not currently working**

### Via Homebrew (for macOS)

Prerequisites:

- [Homebrew](https://brew.sh/)

```
brew install danielrivasmd/Bender
```



<!-- TODO: -->
**Not currently working**

### Via APT (for Debian-based Linux distros)

```
curl -SsL https://fbecart.github.io/ppa/debian/KEY.gpg | sudo apt-key add -
sudo curl -SsL -o /etc/apt/sources.list.d/fbecart.list https://fbecart.github.io/ppa/debian/fbecart.list
sudo apt update
sudo apt install bender
```



<!-- TODO: -->
**Not currently working**

### From Github release

<!-- TODO: -->
**Not currently working**



## Documentation

### Usage

Use `bender -h` or `bender --help` to display help on commandline.

```
"Good news everyone!"
Bender automates Genomic jobs in Slurm systems.
"It's highly addictive!"

Bender creates a convinient command line interphase
with built-in and accessible documentation

Usage:
  bender [command]

Available Commands:
  help                  Help about any command

Flags:
  -h, --help            help for Bender
  -o, --outDir string   Output directory. creates if not exitst

Use "bender [command] --help" for more information about a command.
```



## Examples

<!-- TODO:
add additional example in example folder
 -->

### Example `bender` config

```toml
```
<!-- TODO:
 -->



## _Bender_ for the curious



## Acknowledgements



## License

Bender is distributed under the terms of the GNU GENERAL PUBLIC LICENSE.

See [LICENSE](LICENSE) for details.

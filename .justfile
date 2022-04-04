####################################################################################################

_default:
  @just --list

####################################################################################################

# print justfile
print:
  bat .justfile --language make

####################################################################################################

# edit justfile
edit:
  micro .justfile

####################################################################################################

# aliases
alias b := buildLinux
alias i := install
alias p := hermesPawsey
alias u := hermesUppmax

####################################################################################################

# build bender for OSX & store `excalibur`
buildOSX:
  #!/bin/bash
  set -euo pipefail

  # declarations
  source .just.sh

  echo "Building..."
  go build -v -o ${bender}/excalibur/bender

####################################################################################################

# build bender for linux & store `excalibur`
buildLinux:
  #!/bin/bash
  set -euo pipefail

  # declarations
  source .just.sh

  echo "Building..."
  env GOOS=linux GOARCH=amd64 go build -v -o ${bender}/excalibur/bender

####################################################################################################

# install bender locally
install:
  #!/bin/bash
  set -euo pipefail

  # declarations
  source .just.sh

  echo "Installing..."
  # Bender
  go install
  mv -v "${HOME}/.go/bin/Bender" "${HOME}/.go/bin/bender"

####################################################################################################

# deliver bender binary & shell scripts Uppmax
hermesUppmax:
  #!/bin/bash
  set -euo pipefail

  # declarations
  source .just.sh

  echo "Deploying to Uppmax..."
  rsync -azvhP "${bender}excalibur/bender" "${uppmaxID}:${uppmaxBin}"

  # link sh scripts
  echo "Linking remotely..."
  rsync -azvhPX "${bender}/sh" "${uppmaxID}:${uppmaxBin}/goTools/"

####################################################################################################

# deliver bender binary & shell scripts Pawsey
hermesPawsey:
  #!/bin/bash
  set -euo pipefail

  # declarations
  source .just.sh

  # transfer binary
  echo "Deploying to Pawsey..."
  rsync -azvhP "${bender}/excalibur/bender" "${pawseyID}:${pawseyBin}/"

  # link sh scripts
  echo "Linking remotely..."
  rsync -azvhPX "${bender}/sh" "${pawseyID}:${pawseyBin}/goTools/"

####################################################################################################

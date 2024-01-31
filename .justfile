####################################################################################################

_default:
  @just --list

####################################################################################################

# print justfile
@show:
  bat .justfile --language make

####################################################################################################

# edit justfile
@edit:
  micro .justfile

####################################################################################################

# aliases

####################################################################################################

# home
home := "${HOME}"
remoteHome := "/home/drivas"
bender := "{{home}}/Factorem/Bender"
gobin := "{{home}}/.go/bin"

#################################################################################

# pawsey
pawseyID := "drivas@topaz.pawsey.org.au"
pawseyBin := "{{remoteHome}}/bin"

# uppmax
uppmaxID := "drivas@rackham.uppmax.uu.se"
uppmaxBin := "{{remoteHome}}/bin"

####################################################################################################

# build bender for OSX & store `excalibur`
osx:
  echo "Building..." && go build -v -o {{bender}}/excalibur/bender

####################################################################################################

# build bender for linux & store `excalibur`
linux:
  echo "Building..." && env GOOS=linux GOARCH=amd64 go build -v -o {{bender}}/excalibur/bender

####################################################################################################

# install bender locally
install:
  echo "Installing..." && go install && mv -v {{gobin}}/Bender {{gobin}}/bender

####################################################################################################

# deliver bender binary & shell scripts Uppmax
hermesUppmax: _deployUppmax _linkUppmax

# transfer binary
_deployUppmax:
  echo "Deploying to Uppmax..." && rsync -azvhP {{bender}}/excalibur/bender {{uppmaxID}}:{{uppmaxBin}}

# link sh scripts
_linkUppmax:
  echo "Linking remotely..." && rsync -azvhPX {{bender}}/sh {{uppmaxID}}:{{uppmaxBin}}/goTools/

####################################################################################################

# deliver bender binary & shell scripts Pawsey
hermesPawsey: _deployPawsey _linkPawsey

# transfer binary
_deployPawsey:
  echo "Deploying to Pawsey..." && rsync -azvhP {{bender}}/excalibur/bender {{pawseyID}}:{{pawseyBin}}

# link sh scripts
_linkPawsey:
  echo "Linking remotely..." && rsync -azvhPX {{bender}}/sh {{pawseyID}}:{{pawseyBin}}/goTools/

####################################################################################################
# compose protocols
####################################################################################################

# build & deploy
pawsey: linux && hermesPawsey
uppmax: linux && hermesUppmax

####################################################################################################

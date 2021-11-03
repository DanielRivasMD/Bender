/*
Copyright © 2021 Daniel Rivas <danielrivasmd@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
package cmd

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"strings"

	"github.com/atrox/homedir"
	"github.com/labstack/gommon/color"
	"github.com/spf13/cobra"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

var (
	library    string
	libraryDir string
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// searchCmd represents the search command
var searchCmd = &cobra.Command{
	Use:   "search",
	Short: "Perform similarity search.",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

Perform similarity search on assemblies
using different tools, i.e., blast, diamond.

First, create database.
Next, execute similarity search.
`,

	Example: `
` + chalk.Cyan.Color("bender") + ` assembly Search ` + chalk.Yellow.Color("diamond") + `
  ` + chalk.Green.Color("--configPath") + ` willLeadToConfig/
  ` + chalk.Green.Color("--configFile") + ` foundConfig.toml
  ` + chalk.Green.Color("--assembly") + ` aweSap01.fa
  ` + chalk.Green.Color("--species") + ` awesomeSapiens`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	ValidArgs: []string{"diamond"},
	Args:      cobra.ExactValidArgs(1),

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// execute logic
		assemblySearch(args[0])

	},

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	assemblyCmd.AddCommand(searchCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	searchCmd.Flags().StringVarP(&library, "library", "l", "", "Library to search against")
	searchCmd.Flags().StringVarP(&libraryDir, "libraryDir", "L", "", "Library directory")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func assemblySearch(searchTool string) {

	// find home directory.
	home, errHomedir := homedir.Dir()
	if errHomedir != nil {
		fmt.Println(errHomedir)
		os.Exit(1)
	}

	// trim suffixes
	libraryT := strings.TrimSuffix(library, ".fasta")
	assemblyT := strings.TrimSuffix(assembly, ".fasta.gz")

	// lineBreaks
	lineBreaks()

	// buffers
	var stdout bytes.Buffer
	var stderr bytes.Buffer

	searchType := ""
	switch searchTool {
	case "diamond":
		searchType = "genomeDiamond.sh"
	}

	// shell call
	commd := home + "/bin/goTools/sh/" + searchType
	shCmd := exec.Command(commd, species, assemblyT, inDir, libraryT, libraryDir, outDir)

	// run
	shCmd.Stdout = &stdout
	shCmd.Stderr = &stderr
	_ = shCmd.Run()

	// stdout
	color.Println(color.Cyan(stdout.String(), color.B))

	// stderr
	if stderr.String() != "" {
		color.Println(color.Red(stderr.String(), color.B))
	}

	// lineBreaks
	lineBreaks()
}

////////////////////////////////////////////////////////////////////////////////////////////////////

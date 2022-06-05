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
	"log"
	"os"
	"os/exec"
	"strings"
	"syscall"

	"github.com/spf13/cobra"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var (
	library     string
	libraryDir  string
	frameshit   string
	blockSize   string
	indexChunks string
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// searchCmd represents the search command
var searchCmd = &cobra.Command{
	Use:   "search",
	Short: "Perform similarity search.",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

Perform similarity search on assemblies
using different tools, i.e., ` + chalk.Yellow.Color("blast") + `, ` + chalk.Yellow.Color("diamond") + `.

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
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	assemblyCmd.AddCommand(searchCmd)

	// flags
	searchCmd.Flags().StringVarP(&library, "library", "l", "", "Library to search against")
	searchCmd.Flags().StringVarP(&libraryDir, "libraryDir", "L", "", "Library directory")
	searchCmd.Flags().StringVarP(&frameshit, "frameshit", "f", "15", "diamond blastx frameshit")
	searchCmd.Flags().StringVarP(&blockSize, "blockSize", "b", "6", "diamond blastx blockSize")
	searchCmd.Flags().StringVarP(&indexChunks, "indexChunks", "x", "1", "diamond blastx indexChunks")

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func assemblySearch(searchTool string) {

	// trim suffixes
	libraryT := strings.TrimSuffix(library, ".fasta")

	// TODO: manipulate sequences to identify ORFs & use blastn
	switch searchTool {

	// diamond
	case "diamond":

		// check executable availability
		command, err := exec.LookPath(searchTool)
		if err != nil {
			panic(err)
		}

		// make database from syncytin protein sequence
		if !fileExist(libraryDir + "/" + libraryT + ".dmnd") {
			log.Fatal("Building diamond database...")
			diamondMakedb := []string{
				searchTool, "makedb",
				"--in", libraryDir + "/" + library,
				"--db", libraryDir + "/" + libraryT + ".dmnd",
			}

			err = syscall.Exec(command, diamondMakedb, os.Environ())
			if err != nil {
				panic(err)
			}
		}

		// use assembly as query with six reading frames
		log.Fatal("Searching diamond database...")
		diamondBlastx := []string{
			searchTool, "blastx",
			"--db", libraryDir + "/" + libraryT + ".dmnd",
			"--query", inDir + "/" + assembly,
			"--frameshift", frameshit,
			"--block-size", blockSize,
			"--index-chunks", indexChunks,
			"--out", outDir + "/" + species + ".tsv",
		}

		err = syscall.Exec(command, diamondBlastx, os.Environ())
		if err != nil {
			panic(err)
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
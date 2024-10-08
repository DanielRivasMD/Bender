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
	"fmt"
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
	frameshit   string
	blockSize   string
	indexChunks string
)

var (
	assemblyT string
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// searchCmd represents the search command
var searchCmd = &cobra.Command{
	Use:   "search",
	Short: "Perform similarity search.",
	Long: chalk.Green.Color("Daniel Rivas <danielrivasmd@gmail.com>") + `

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

	Run: func(κ *cobra.Command, args []string) {

		// execute logic
		assemblySearch(args[0])

	},
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	assemblyCmd.AddCommand(searchCmd)

	// flags
	searchCmd.Flags().StringVarP(&frameshit, "frameshit", "f", "15", "diamond blastx frameshit")
	searchCmd.Flags().StringVarP(&blockSize, "blockSize", "b", "2", "diamond blastx blockSize")
	searchCmd.Flags().StringVarP(&indexChunks, "indexChunks", "x", "1", "diamond blastx indexChunks")

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// search assemblies with external tool
func assemblySearch(searchTool string) {

	// trim suffixes
	libraryT := strings.TrimSuffix(library, ".fasta")

	// accommodate different fasta formats
	switch {
	case strings.Contains(assembly, "fna"):
		assemblyT = strings.TrimSuffix(assembly, ".fna.gz")
	case strings.Contains(assembly, "fasta"):
		assemblyT = strings.TrimSuffix(assembly, ".fasta.gz")
	default:
		log.Fatal("fasta format not recognized")
	}

	// trim prefix
	scaffoldT := strings.TrimPrefix(scaffold, assemblyT+"_")

	// use search tool
	switch searchTool {

	// diamond
	case "diamond":

		// check executable availability
		κ, ε := exec.LookPath(searchTool)
		if ε != nil {
			panic(ε)
		}

		// use assembly as query with six reading frames
		fmt.Println("Searching diamond database...")
		diamondBlastx := []string{
			searchTool, "blastx",
			"--db", libraryDir + "/" + libraryT + ".dmnd",
			"--query", inDir + "/" + assemblyT + "_" + scaffoldT,
			"--frameshift", frameshit,
			"--block-size", blockSize,
			"--index-chunks", indexChunks,
			"--out", outDir + "/" + species + "_" + scaffoldT + ".tsv",
		}

		// execute command
		ε = syscall.Exec(κ, diamondBlastx, os.Environ())
		if ε != nil {
			panic(ε)
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////

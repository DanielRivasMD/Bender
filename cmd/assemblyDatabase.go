/*
Copyright © 2022 Daniel Rivas <danielrivasmd@gmail.com>

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
	"os"
	"os/exec"
	"strings"
	"syscall"

	"github.com/spf13/cobra"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var ()

////////////////////////////////////////////////////////////////////////////////////////////////////

// databaseCmd represents the database command
var databaseCmd = &cobra.Command{
	Use:   "database",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	ValidArgs: []string{"diamond"},
	Args:      cobra.ExactValidArgs(1),

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(κ *cobra.Command, args []string) {

		// execute logic
		assemblyDatabase(args[0])

	},
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	assemblyCmd.AddCommand(databaseCmd)

	// flags

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func assemblyDatabase(searchTool string) {

	// trim suffixes
	libraryT := strings.TrimSuffix(library, ".fasta")

	// TODO: manipulate sequences to identify ORFs & use blastn
	switch searchTool {

	// diamond
	case "diamond":

		// check executable availability
		κ, ε := exec.LookPath(searchTool)
		if ε != nil {
			panic(ε)
		}

		// make database from syncytin protein sequence
		if !fileExist(libraryDir + "/" + libraryT + ".dmnd") {
			fmt.Println("Building diamond database...")
			diamondMakedb := []string{
				searchTool, "makedb",
				"--in", libraryDir + "/" + library,
				"--db", libraryDir + "/" + libraryT + ".dmnd",
			}

			ε = syscall.Exec(κ, diamondMakedb, os.Environ())
			if ε != nil {
				panic(ε)
			}
		}

	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

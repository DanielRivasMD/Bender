/*
Copyright © 2020 Daniel Rivas <danielrivasmd@gmail.com>

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
	"log"
	"os"
	"os/exec"

	"github.com/atrox/homedir"
	"github.com/labstack/gommon/color"
	"github.com/spf13/cobra"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

var (
	assemblyBlast    string
	assemblyDirBlast string
	speciesBlast     string
	libraryBlast     string
	libraryDirBlast  string
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// GenomeBlastCmd represents the GenomeBlast command
var GenomeBlastCmd = &cobra.Command{
	Use:   "GenomeBlast",
	Short: "Blast an assembly with customized sequence library",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

GenomeBlast performs several operations:
- Inputs a customized sequence
- Creates a database from an assembly
- Searches possible homology in an assembly
- Formats output the values`,
	Example: `
bender GenomeBlast -s SuperMouse -a toSearch.fa.gz -A findSearch -l toQuery.fa -L findQuery`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// Find home directory.
		home, errHomedir := homedir.Dir()
		if errHomedir != nil {
			fmt.Println(errHomedir)
			os.Exit(1)
		}

		// lineBreaks
		lineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

		// shell call
		commd := home + "/bin/goTools/sh/GenomeBlast.sh"
		shCmd := exec.Command(commd, speciesBlast, assemblyBlast, assemblyDirBlast, libraryBlast, libraryDirBlast, outDir)

		// run
		shCmd.Stdout = &stdout
		shCmd.Stderr = &stderr
		errShCmd := shCmd.Run()

		if errShCmd != nil {
			log.Printf("error: %v\n", errShCmd)
		}

		// stdout
		color.Println(color.Cyan(stdout.String(), color.B))

		// stderr
		if stderr.String() != "" {
			color.Println(color.Red(stderr.String(), color.B))
		}

		blastFilter(outDir + speciesBlast)

		// lineBreaks
		lineBreaks()

	},

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	rootCmd.AddCommand(GenomeBlastCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	GenomeBlastCmd.Flags().StringVarP(&speciesBlast, "species", "s", "", "Species")
	GenomeBlastCmd.Flags().StringVarP(&assemblyBlast, "assembly", "a", "", "Assembly to BLAST")
	GenomeBlastCmd.Flags().StringVarP(&assemblyDirBlast, "assemblyDir", "A", "", "Assembly directory")
	GenomeBlastCmd.Flags().StringVarP(&libraryBlast, "library", "l", "", "Library to BLAST against")
	GenomeBlastCmd.Flags().StringVarP(&libraryDirBlast, "libraryDir", "L", "", "Library directory")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

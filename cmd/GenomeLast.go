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

	"github.com/atrox/homedir"
	"github.com/labstack/gommon/color"
	"github.com/spf13/cobra"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

var (
	assemblyLast    string
	assemblyDirLast string
	speciesLast     string
	libraryLast     string
	libraryDirLast  string
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// GenomeLastCmd represents the GenomeLast command
var GenomeLastCmd = &cobra.Command{
	Use:   "GenomeLast",
	Short: "Last genetic sequences",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

GenomeLast performs...
TODO: fill up
`,

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
		commd := home + "/bin/goTools/sh/GenomeLast.sh"
		shCmd := exec.Command(commd, speciesLast, assemblyLast, assemblyDirLast, libraryLast, libraryDirLast, outDir)

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

		// TODO: filter records after check
		// // filter records

		// lineBreaks
		lineBreaks()

	},

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

func init() {
	rootCmd.AddCommand(GenomeLastCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	GenomeLastCmd.Flags().StringVarP(&speciesLast, "species", "s", "", "Species")
	GenomeLastCmd.Flags().StringVarP(&assemblyLast, "assembly", "a", "", "Assembly to Last")
	GenomeLastCmd.Flags().StringVarP(&assemblyDirLast, "assemblyDir", "A", "", "Assembly directory")
	GenomeLastCmd.Flags().StringVarP(&libraryLast, "library", "l", "", "Library to Last against")
	GenomeLastCmd.Flags().StringVarP(&libraryDirLast, "libraryDir", "L", "", "Library directory")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

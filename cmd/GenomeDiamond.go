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
	assemblyDiamond    string
	assemblyDirDiamond string
	speciesDiamond     string
	libraryDiamond     string
	libraryDirDiamond  string
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// GenomeDiamondCmd represents the GenomeDiamond command
var GenomeDiamondCmd = &cobra.Command{
	Use:   "GenomeDiamond",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,

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
		commd := home + "/bin/goTools/sh/GenomeDiamond.sh"
		shCmd := exec.Command(commd, speciesDiamond, assemblyDiamond, assemblyDirDiamond, libraryDiamond, libraryDirDiamond, outDir)

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

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	rootCmd.AddCommand(GenomeDiamondCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	GenomeDiamondCmd.Flags().StringVarP(&speciesDiamond, "species", "s", "", "Species")
	GenomeDiamondCmd.Flags().StringVarP(&assemblyDiamond, "assembly", "a", "", "Assembly to Diamond")
	GenomeDiamondCmd.Flags().StringVarP(&assemblyDirDiamond, "assemblyDir", "A", "", "Assembly directory")
	GenomeDiamondCmd.Flags().StringVarP(&libraryDiamond, "library", "l", "", "Library to Diamond against")
	GenomeDiamondCmd.Flags().StringVarP(&libraryDirDiamond, "libraryDir", "L", "", "Library directory")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

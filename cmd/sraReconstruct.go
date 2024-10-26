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
	"log"
	"os/exec"

	"github.com/labstack/gommon/color"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var ()

////////////////////////////////////////////////////////////////////////////////////////////////////

// reconstructCmd represents the reconstruct command
var reconstructCmd = &cobra.Command{
	Use:   "reconstruct",
	Short: "Reconstruct SRA binaries.",
	Long: chalk.Green.Color("Daniel Rivas <danielrivasmd@gmail.com>") + `

Reconstruct binary SRA files to fastq.
`,

	Example: `
` + chalk.Cyan.Color("bender") + ` sra reconstruct --inDir projPath/ --outDir ReconstructPath --file coolBinaryList.txt --split-files`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(κ *cobra.Command, args []string) {

		// flags
		inDir, _ := κ.Flags().GetString("inDir")
		outDir, _ := κ.Flags().GetString("outDir")
		γ, _ := κ.Flags().GetString("verbose")
		ƒ, _ := κ.Flags().GetString("file")

		splitFiles, _ := κ.Flags().GetString("split-files")

		// lineBreaks
		lineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

		// shell call
		commd := findHome() + "/bin/goTools/sh/SRAreconstruct.sh"
		shCmd := exec.Command(commd, inDir, outDir, γ, ƒ, splitFiles)

		// run
		shCmd.Stdout = &stdout
		shCmd.Stderr = &stderr
		err := shCmd.Run()

		if err != nil {
			log.Printf("error: %v\n", err)
		}

		// stdout
		color.Println(color.Cyan(stdout.String(), color.B))

		// stderr
		if stderr.String() != "" {
			color.Println(color.Red(stderr.String(), color.B))
		}

		// lineBreaks
		lineBreaks()

	},
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	sraCmd.AddCommand(reconstructCmd)

	// TODO: update flag bindings
	// flags
	reconstructCmd.Flags().StringP("file", "f", "", "File containing binary list")
	reconstructCmd.MarkFlagRequired("file")
	viper.BindPFlag("file", reconstructCmd.Flags().Lookup("file"))

	// TODO: toggle?
	reconstructCmd.Flags().StringP("split-files", "s", "false", "Indicates whether binary will be reconstructed in different files")
	reconstructCmd.MarkFlagRequired("split-files")
	viper.BindPFlag("split-files", reconstructCmd.Flags().Lookup("split-files"))

}

////////////////////////////////////////////////////////////////////////////////////////////////////

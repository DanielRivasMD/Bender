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
	"log"
	"os"
	"os/exec"

	"github.com/DanielRivasMD/Bender/aux"
	"github.com/atrox/homedir"
	"github.com/labstack/gommon/color"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

// ReconstructCmd represents the Reconstruct command
var ReconstructCmd = &cobra.Command{
	Use:   "Reconstruct",
	Short: "Reconstruct SRA binaries",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

SRA reconstruct binary SRA files to fastq
`,
	Example: `
Bender SRA Reconstruct --inDir projPath/ --outDir ReconstructPath --file coolBinaryList.txt --split-files`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// find home directory.
		home, errHomedir := homedir.Dir()
		if errHomedir != nil {
			fmt.Println(errHomedir)
			os.Exit(1)
		}

		// flags
		inDir, _ := cmd.Flags().GetString("inDir")
		outDir, _ := cmd.Flags().GetString("outDir")
		verbose, _ := cmd.Flags().GetString("verbose")
		file, _ := cmd.Flags().GetString("file")

		splitFiles, _ := cmd.Flags().GetString("split-files")

		// lineBreaks
		aux.LineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

			// shell call
			commd := home + "/bin/goTools/sh/SRAreconstruct.sh"
			shCmd := exec.Command(commd, inDir, outDir, verbose, file, splitFiles)

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
		aux.LineBreaks()

	},

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

func init() {
	SRACmd.AddCommand(ReconstructCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	ReconstructCmd.Flags().StringP("file", "f", "", "File containing binary list")
	ReconstructCmd.MarkFlagRequired("file")
	viper.BindPFlag("file", ReconstructCmd.Flags().Lookup("file"))

	// TODO: toggle?
	ReconstructCmd.Flags().StringP("split-files", "s", "false", "Indicates whether binary will be reconstructed in different files")
	ReconstructCmd.MarkFlagRequired("split-files")
	viper.BindPFlag("split-files", ReconstructCmd.Flags().Lookup("split-files"))

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

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
	"os/exec"
	"strings"

	"github.com/labstack/gommon/color"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var ()

////////////////////////////////////////////////////////////////////////////////////////////////////

// fastqExtractCmd represents the bamFastqExtract command
var fastqExtractCmd = &cobra.Command{
	Use:   "fastqExtract",
	Short: "Extract FASTQ from BAM files.",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

Dissect BAM files and retrieve FASTQ files.
`,

	Example: `
` + chalk.Cyan.Color("bender") + ` bamFastqExtract -f aCloneOfMyOwn.bam`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(κ *cobra.Command, args []string) {

		// flags
		storageDir, _ := κ.Flags().GetString("outDir")

		ƒ, _ := κ.Flags().GetString("file")
		ƒ = strings.TrimSuffix(ƒ, ".bam")

		directory, _ := κ.Flags().GetString("inDir")

		ɣ, _ := κ.Flags().GetString("ɣ")

		// lineBreaks
		lineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

		// shell call
		commd := findHome() + "/bin/goTools/sh/bamFastqExtract.sh"
		shCmd := exec.Command(commd, ƒ, directory, ɣ, storageDir)

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

	},
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	bamCmd.AddCommand(fastqExtractCmd)

	// flags
	fastqExtractCmd.Flags().StringP("file", "f", "", "Alignment file. bam format")
	fastqExtractCmd.MarkFlagRequired("file")
	viper.BindPFlag("file", fastqExtractCmd.Flags().Lookup("file"))

}

////////////////////////////////////////////////////////////////////////////////////////////////////

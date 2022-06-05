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
	"strings"

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
		ƒ, _ := κ.Flags().GetString("file")
		ƒ = strings.TrimSuffix(ƒ, ".bam")

		directory, _ := κ.Flags().GetString("inDir")

		ɣ, _ := κ.Flags().GetString("verbose")

		outDir, _ := κ.Flags().GetString("outDir")

		// execute shell
		execCmd("bamFastqExtract.sh", ƒ, directory, ɣ, outDir)

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

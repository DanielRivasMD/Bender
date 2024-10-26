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
	"strings"

	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var ()

////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO: update bwa cmd
// bwaCmd represents the bwa command
var bwaCmd = &cobra.Command{
	Use:   "bwa",
	Short: "Align fasta using BWA.",
	Long: chalk.Green.Color("Daniel Rivas <danielrivasmd@gmail.com>") + `

` + chalk.Green.Color("Bender") + ` will align fasta files to a reference genome.
Additionally, ` + chalk.Green.Color("Bender") + ` perform quality control prealignment through FastX.
`,

	Example: `
` + chalk.Cyan.Color("bender") + ` bwa -f lastOneAlive.fa -r cyclops.fa`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(κ *cobra.Command, args []string) {

		// flags
		ƒ, _ := κ.Flags().GetString("file")
		ƒ = strings.TrimSuffix(ƒ, ".fa")

		directory, _ := κ.Flags().GetString("inDir")

		ξ, _ := κ.Flags().GetString("verbose")

		outDir, _ := κ.Flags().GetString("outDir")

		// execute shell
		execCmd("bwa.sh", ƒ, directory, ξ, outDir)

	},
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	rootCmd.AddCommand(bwaCmd)

	// flags
	bwaCmd.Flags().StringP("file", "f", "", "File to align. fasta format")
	bwaCmd.MarkFlagRequired("file")
	viper.BindPFlag("file", bwaCmd.Flags().Lookup("file"))

	bwaCmd.Flags().StringP("reference", "r", "", "Reference file. fasta format")
	bwaCmd.MarkFlagRequired("reference")
	viper.BindPFlag("reference", bwaCmd.Flags().Lookup("reference"))

}

////////////////////////////////////////////////////////////////////////////////////////////////////

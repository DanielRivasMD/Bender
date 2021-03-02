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
	"strings"

	"github.com/DanielRivasMD/Bender/auxiliary"
	"github.com/atrox/homedir"
	"github.com/labstack/gommon/color"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

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

		// bound flags
		outDir := viper.GetString("outDir")

		genome := viper.GetString("genome")
		if genome == "" {
			genome, _ = cmd.Flags().GetString("genome")
		}
		genome = strings.TrimSuffix(genome, ".fasta")

		genomeDir := viper.GetString("genomeDir")

		library := viper.GetString("library")

		libraryDir := viper.GetString("libraryDir")

		// lineBreaks
		auxiliary.LineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

		// shell call
		commd := home + "/bin/goTools/sh/GenomeLast.sh"
		shCmd := exec.Command(commd, genome, genomeDir, library, libraryDir, outDir)

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

		// auxiliary.FileReader(outDir + genome)

		// lineBreaks
		auxiliary.LineBreaks()

	},

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

func init() {
	rootCmd.AddCommand(GenomeLastCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	GenomeLastCmd.Flags().StringP("genome", "g", "", "Genome to LAST")
	viper.BindPFlag("genome", GenomeLastCmd.Flags().Lookup("genome"))

	GenomeLastCmd.Flags().StringP("genomeDir", "G", "", "Genome directory")
	viper.BindPFlag("genomeDir", GenomeLastCmd.Flags().Lookup("genomeDir"))

	GenomeLastCmd.Flags().StringP("library", "l", "", "Library to LAST against")
	viper.BindPFlag("library", GenomeLastCmd.Flags().Lookup("library"))

	GenomeLastCmd.Flags().StringP("libraryDir", "L", "", "Library directory")
	viper.BindPFlag("libraryDir", GenomeLastCmd.Flags().Lookup("libraryDir"))

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

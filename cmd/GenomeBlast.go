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
	"strings"

	"github.com/DanielRivasMD/Bender/aux"
	"github.com/atrox/homedir"
	"github.com/labstack/gommon/color"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

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
Bender GenomeBlast -l toQuery.fa -L findQuery/ -g toSearch.fa -G findSearch/`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// Find home directory.
		home, errHomedir := homedir.Dir()
		if errHomedir != nil {
			fmt.Println(errHomedir)
			os.Exit(1)
		}

		// flags
		storageDir, _ := cmd.Flags().GetString("outDir")

		genome, _ := cmd.Flags().GetString("genome")
		genome = strings.TrimSuffix(genome, ".fasta")

		genomeDir, _ := cmd.Flags().GetString("genomeDir")

		library, _ := cmd.Flags().GetString("library")

		libraryDir, _ := cmd.Flags().GetString("libraryDir")

		// lineBreaks
		aux.LineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

		// shell call
		commd := home + "/bin/goTools/sh/GenomeBlast.sh"
		shCmd := exec.Command(commd, genome, genomeDir, library, libraryDir, storageDir)

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

		aux.FileReader(storageDir + genome)

		// lineBreaks
		aux.LineBreaks()

	},

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

func init() {
	rootCmd.AddCommand(GenomeBlastCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	GenomeBlastCmd.Flags().StringP("genome", "g", "", "Genome to BLAST")
	GenomeBlastCmd.MarkFlagRequired("genome")
	viper.BindPFlag("genome", GenomeBlastCmd.Flags().Lookup("genome"))

	GenomeBlastCmd.Flags().StringP("genomeDir", "G", "", "Genome directory")
	GenomeBlastCmd.MarkFlagRequired("genomeDir")
	viper.BindPFlag("genomeDir", GenomeBlastCmd.Flags().Lookup("genomeDir"))

	GenomeBlastCmd.Flags().StringP("library", "l", "", "Library to BLAST against")
	GenomeBlastCmd.MarkFlagRequired("library")
	viper.BindPFlag("library", GenomeBlastCmd.Flags().Lookup("library"))

	GenomeBlastCmd.Flags().StringP("libraryDir", "L", "", "Library directory")
	GenomeBlastCmd.MarkFlagRequired("libraryDir")
	viper.BindPFlag("libraryDir", GenomeBlastCmd.Flags().Lookup("libraryDir"))

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

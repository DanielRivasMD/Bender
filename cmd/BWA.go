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

	"github.com/atrox/homedir"
	"github.com/labstack/gommon/color"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

// BWACmd represents the BWA command
var BWACmd = &cobra.Command{
	Use:   "BWA",
	Short: "Align fasta using BWA",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

BWA aligns fasta files to a reference genome.
Additionally, BWA can perform quality control
prealignment through FastX`,
	Example: `
Bender BWA -f lastOneAlive.fa -r cyclops.fa`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// find home directory.
		home, errHomedir := homedir.Dir()
		if errHomedir != nil {
			fmt.Println(errHomedir)
			os.Exit(1)
		}

		// flags
		storageDir, _ := cmd.Flags().GetString("outDir")

		file, _ := cmd.Flags().GetString("file")
		file = strings.TrimSuffix(file, ".fa")

		directory, _ := cmd.Flags().GetString("inDir")

		verbose, _ := cmd.Flags().GetString("verbose")

		// lineBreaks
		lineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

		// shell call
		commd := home + "/bin/goTools/sh/BWA.sh"
		shCmd := exec.Command(commd, file, directory, verbose, storageDir)

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

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

func init() {
	rootCmd.AddCommand(BWACmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	BWACmd.Flags().StringP("file", "f", "", "File to align. fasta format")
	BWACmd.MarkFlagRequired("file")
	viper.BindPFlag("file", BWACmd.Flags().Lookup("file"))

	BWACmd.Flags().StringP("reference", "r", "", "Reference file. fasta format")
	BWACmd.MarkFlagRequired("reference")
	viper.BindPFlag("reference", BWACmd.Flags().Lookup("reference"))

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

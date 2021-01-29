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

	"github.com/DanielRivasMD/Bender/aux"
	"github.com/atrox/homedir"
	"github.com/labstack/gommon/color"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

// HTSeqCmd represents the HTSeq command
var HTSeqCmd = &cobra.Command{
	Use:   "HTSeq",
	Short: "Wrapper for HTSeq",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

HTSeq wraps the namesake python package,
which is a command line tool application for
processing high-throughput sequencing data`,
	Example: `
Bender HTSeq -f forDESeq`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// find home directory.
		home, errHomedir := homedir.Dir()
		if errHomedir != nil {
			fmt.Println(errHomedir)
			os.Exit(1)
		}

		// flags
		file, _ := cmd.Flags().GetString("file")

		verbose, _ := cmd.Flags().GetString("verbose")

		// lineBreaks
		aux.LineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

		// shell call
		commd := home + "/bin/goTools/sh/HTSeq.sh"
		shCmd := exec.Command(commd, file, verbose)

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
	rootCmd.AddCommand(HTSeqCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	HTSeqCmd.Flags().StringP("file", "f", "", "File")
	HTSeqCmd.MarkFlagRequired("file")
	viper.BindPFlag("file", HTSeqCmd.Flags().Lookup("file"))

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

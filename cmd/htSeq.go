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

	"github.com/labstack/gommon/color"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var ()

////////////////////////////////////////////////////////////////////////////////////////////////////

// htSeqCmd represents the htseq command
var htSeqCmd = &cobra.Command{
	Use:   "htSeq",
	Short: "Wrapper for HTseq.",
	Long: chalk.Green.Color("Daniel Rivas <danielrivasmd@gmail.com>") + `

Wrap HTseq python package,
a command line tool application for
processing high-throughput sequencing data.
`,

	Example: `
` + chalk.Cyan.Color("bender") + ` htseq -f forDESeq`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(κ *cobra.Command, args []string) {

		// flags
		ƒ, _ := κ.Flags().GetString("file")

		ɣ, _ := κ.Flags().GetString("verbose")

		// lineBreaks
		lineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

		// shell call
		commd := findHome() + "/bin/goTools/sh/htSeq.sh"
		shCmd := exec.Command(commd, ƒ, ɣ)

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
	rootCmd.AddCommand(htSeqCmd)

	// flags
	htSeqCmd.Flags().StringP("file", "f", "", "File")
	htSeqCmd.MarkFlagRequired("file")
	viper.BindPFlag("file", htSeqCmd.Flags().Lookup("file"))

}

////////////////////////////////////////////////////////////////////////////////////////////////////

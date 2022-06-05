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

// repeatModelerCmd represents the repeatModeler command
var repeatModelerCmd = &cobra.Command{
	Use:   "repeatModeler",
	Short: "Create Repeat Library from assembly.",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

Create Repeat Library from assembly using RepeatModeler v2.0.1.

First, build a database.
Next, create a libray.
`,

	Example: `
` + chalk.Cyan.Color("bender") + ` repeatModeler -r phillipFry.fa`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(κ *cobra.Command, args []string) {

		// flags
		storageDir, _ := κ.Flags().GetString("outDir")

		directory, _ := κ.Flags().GetString("inDir")

		ɣ, _ := κ.Flags().GetString("verbose")

		// bound flags
		reference := viper.GetString("reference")
		reference = strings.TrimSuffix(reference, ".fa")

		// lineBreaks
		lineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

		// shell call
		commd := findHome() + "/bin/goTools/sh/repeatModeler.sh"
		shCmd := exec.Command(commd, reference, directory, ɣ, storageDir)

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
	rootCmd.AddCommand(repeatModelerCmd)

	// flags
	repeatModelerCmd.Flags().StringP("reference", "r", "", "Reference file. fasta format")
	viper.BindPFlag("reference", repeatModelerCmd.Flags().Lookup("reference"))

}

////////////////////////////////////////////////////////////////////////////////////////////////////

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
	"github.com/ttacon/chalk"
)

// fetchCmd represents the fetch command
var fetchCmd = &cobra.Command{
	Use:   "fetch",
	Short: "Fetch SRA accessions.",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

Collect files from Short Read Archive (SRA)
and check the state of the downloads.
`,

	Example: `
` + chalk.Cyan.Color("bender") + ` sra fetch --inDir projPath/ --outDir logPath --file coolFileList.txt --iterations 20 --max-size 35G`,

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

		file = strings.TrimSuffix(file, ".txt")
		maxIt, _ := cmd.Flags().GetString("iterations")
		maxSize, _ := cmd.Flags().GetString("max-size")

		// lineBreaks
		lineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

		// shell call
		commd := home + "/bin/goTools/sh/fetch.sh"
		shCmd := exec.Command(commd, inDir, outDir, verbose, file, maxIt, maxSize, verbose)

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
	sraCmd.AddCommand(fetchCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	fetchCmd.Flags().StringP("file", "f", "", "File containing accession numbers")
	fetchCmd.MarkFlagRequired("file")
	viper.BindPFlag("file", fetchCmd.Flags().Lookup("file"))

	fetchCmd.Flags().StringP("max-size", "m", "20G", "File size limit to download")
	fetchCmd.MarkFlagRequired("max-size")
	viper.BindPFlag("max-size", fetchCmd.Flags().Lookup("max-size"))

	fetchCmd.Flags().StringP("iterations", "t", "10", "Number of iterations for download to fail")
	fetchCmd.MarkFlagRequired("iterations")
	viper.BindPFlag("iterations", fetchCmd.Flags().Lookup("iterations"))

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

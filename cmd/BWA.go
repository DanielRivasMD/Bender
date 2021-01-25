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

// BWACmd represents the BWA command
var BWACmd = &cobra.Command{
	Use:   "BWA",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// Find home directory.
		home, errHomedir := homedir.Dir()
		if errHomedir != nil {
			fmt.Println(errHomedir)
			os.Exit(1)
		}

		// read config
		var config Config
		errConf := viper.Unmarshal(&config)
		if errConf != nil {
			log.Fatalf("could not decode config into struct: %v", errConf)
		}

		// Flags
		storageDir := config.StorageDir

		file, _ := cmd.Flags().GetString("file")
		file = strings.TrimSuffix(file, ".fa")

		directory, _ := cmd.Flags().GetString("directory")

		verbose, _ := cmd.Flags().GetString("verbose")

		// lineBreaks
		aux.LineBreaks()

		// buffers
		var stdout bytes.Buffer
		var stderr bytes.Buffer

		// shell call
		commd := home + "/bin/goTools/sh/BWA.sh"
		shCmd := exec.Command(commd, file, directory, verbose, storageDir)

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
	rootCmd.AddCommand(BWACmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// Flags
	BWACmd.Flags().StringP("file", "f", "", "file to align. fasta format")
	BWACmd.MarkFlagRequired("file")

	BWACmd.Flags().StringP("directory", "d", "", "directory")
	BWACmd.MarkFlagRequired("directory")

	BWACmd.Flags().StringP("verbose", "v", "false", "verbosity")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

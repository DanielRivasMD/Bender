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
	"fmt"
	"log"
	"os"

	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var cfgFile string

////////////////////////////////////////////////////////////////////////////////////////////////////

var rootCmd = &cobra.Command{
	Use:   "Bender",
	Short: "A robot to handle Slurm Genomic jobs",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

"Good news everyone!"
Bender is a robot for automation on
Genomic jobs in Slurm systems.
"It's highly addictive!"

Bender creates a convinient command line interphase
with built-in and accessible documentation`,
	Version: "v0.1",
	Example: `
Bender help`,
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func Execute() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	cobra.OnInitialize(initConfig)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// persistent flags
	rootCmd.PersistentFlags().StringVarP(&cfgFile, "config", "c", "", "Config file")

	rootCmd.PersistentFlags().StringP("inDir", "i", "", "Directory where input files are located")
	rootCmd.MarkFlagRequired("inDir")
	viper.BindPFlag("inDir", rootCmd.PersistentFlags().Lookup("inDir"))

	rootCmd.PersistentFlags().StringP("outDir", "o", "", "Output directory. Creates if not exitst")
	rootCmd.MarkFlagRequired("outDir")
	viper.BindPFlag("outDir", rootCmd.PersistentFlags().Lookup("outDir"))

	rootCmd.PersistentFlags().StringP("verbose", "v", "false", "Verbosity switch")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func initConfig() {

	////////////////////////////////////////////////////////////////////////////////////////////////////

	if cfgFile != "" {
		// use config file from the flag.
		viper.SetConfigFile(cfgFile)

		// error handler
		err := viper.ReadInConfig()
		if err != nil {
			log.Fatalf("could not read config file: %v", err)
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

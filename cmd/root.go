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

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "Bender",
	Short: "A robot to handle Slurm Genomic jobs",
	Long: `"Good news everyone!"
Bender automates Genomic jobs in Slurm systems.
"It's highly addictive!"

Bender creates a convinient command line interphase
with built-in and accessible documentation`,
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}

func init() {
	cobra.OnInitialize(initConfig)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// Flags
	rootCmd.PersistentFlags().StringP("outDir", "o", "", "Output directory. creates if not exitst")
	viper.BindPFlag("outDir", rootCmd.PersistentFlags().Lookup("outDir"))

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// initialize configuration
var (
	defaults = map[string]interface{}{
		"outDir":      "data/output/",
		"assemblyDir": "data/DNAzoo/",
		"libraryDir":  "data/syncytinDB/",
	}

	// extension is autodetected
	configName = "bender"

	// search only at local directory
	configPaths = []string{
		".",
	}
)

// load config on struct
type Config struct {
	OutDir      string
	AssemblyDir string
	LibraryDir  string
}

func initConfig() {

	// set defaults
	for k, v := range defaults {
		viper.SetDefault(k, v)
	}

	// set config name
	viper.SetConfigName(configName)
	for _, p := range configPaths {
		viper.AddConfigPath(p)
	}

	// error handler
	err := viper.ReadInConfig()
	if err != nil {
		log.Fatalf("could not read config file: %v", err)
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////

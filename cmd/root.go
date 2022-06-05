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
	"os"
	"strings"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
	"github.com/spf13/viper"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var (
	cfgPath     string
	cfgFile     string
	verboseBool string
	inDir       string
	outDir      string
)

////////////////////////////////////////////////////////////////////////////////////////////////////

var rootCmd = &cobra.Command{
	Use:     "bender",
	Version: "v0.3",
	Short:   "A robot to handle Slurm Genomic jobs.",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

"Good news everyone!"
` + chalk.Green.Color("Bender") + chalk.Magenta.Color(` is a robot for automation on
Genomic jobs in Slurm systems.
`) + chalk.Yellow.Color("It's highly addictive!") + `

` + chalk.Green.Color("Bender") + ` creates a convenient command line interphase
with built-in and accessible documentation.
`,

	Example: `
` + chalk.Cyan.Color("bender") + ` help`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	PersistentPreRun: func(κ *cobra.Command, arg []string) {

		initializeConfig(κ, cfgPath, strings.TrimSuffix(cfgFile, ".toml"))

	},
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func Execute() {
	ε := rootCmd.Execute()
	if ε != nil {
		log.Fatal(ε)
		os.Exit(1)
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func initializeConfig(κ *cobra.Command, configPath string, configName string) error {

	// initialize viper
	ʌ := viper.New()

	// collect config path & file from persistent flags
	ʌ.AddConfigPath(configPath)
	ʌ.SetConfigName(configName)

	// read the config file
	ε := ʌ.ReadInConfig()
	if ε != nil {
		// okay if there isn't a config file
		_, ϙ := ε.(viper.ConfigFileNotFoundError)
		if !ϙ {
			// return an error if we cannot parse the config file
			return ε
		}
	}

	// bind flags to viper
	bindFlags(κ, ʌ)

	return nil
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// bind each cobra flag to its associated viper configuration
func bindFlags(κ *cobra.Command, ʌ *viper.Viper) {

	κ.Flags().VisitAll(func(σ *pflag.Flag) {

		// apply the viper config value to the flag when the flag is not set and viper has a value
		if !σ.Changed && ʌ.IsSet(σ.Name) {
			ν := ʌ.Get(σ.Name)
			κ.Flags().Set(σ.Name, fmt.Sprintf("%v", ν))
		}
	})
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// fileExist checks if a file exists and is not a directory before try using it to prevent further errors
func fileExist(ƒ string) bool {
	info, ε := os.Stat(ƒ)
	if os.IsNotExist(ε) {
		return false
	}
	return !info.IsDir()
}

////////////////////////////////////////////////////////////////////////////////////////////////////
func init() {

	// persistent flags
	rootCmd.PersistentFlags().StringVarP(&cfgFile, "configFile", "c", "", "Config file")
	rootCmd.PersistentFlags().StringVarP(&cfgPath, "configPath", "C", ".", "Path to config file")
	rootCmd.PersistentFlags().StringVarP(&inDir, "inDir", "I", ".", "Directory where input files are located")
	rootCmd.PersistentFlags().StringVarP(&outDir, "outDir", "O", ".", "Output directory. Creates if not exitst")
	rootCmd.PersistentFlags().StringVarP(&verboseBool, "verbose", "v", "false", "Verbosity switch")

}

////////////////////////////////////////////////////////////////////////////////////////////////////

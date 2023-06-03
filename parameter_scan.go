package main

import (
	"fmt"
	"os"
	"flag"
	"errors"
	"strings"
	"sync"
	fp "path/filepath"
	ex "os/exec"
	"regexp"
)

func check(e error) {
    if e != nil {
        panic(e)
    }
}

func create_dir(d string) {
	if _, err := os.Stat(d); errors.Is(err, os.ErrNotExist) {
		os.Mkdir(d, 0700) //os.Mkdir(d, os.ModeDir)
	}
}

func run_hpair(in string, out string) {
    exe := "./run"

	fmt.Println("Command: " + exe + " " + in + " " + out)
    cmd := ex.Command(exe, in, out)

    _, err := cmd.Output()
    if err != nil {
        fmt.Println(err.Error())
        return
    }
}

func create_input_file(tri string, template []byte, out string) {
	var tri2 string;
	if strings.Contains(tri, "m") {
		tri2 = strings.Replace(tri, "m", "-", -1)
	} else {
		tri2 = tri
	}
	
	template_lines := strings.Split(string(template), "\n")
	for i, line := range template_lines {
		if strings.Contains(line, "C_3") {
			template_lines[i] = "C_3      =  " + tri2 + ".0D0";
		}
	}

	output_lines := strings.Join(template_lines, "\n")
	err := os.WriteFile(out, []byte(output_lines), 0644)
	check(err)
}

func extract_cross_sections(out string) ([]string, []string) {
	f, err := os.ReadFile(out)
	check(err)

	var matches_born [][]string;
	var matches_nlo [][]string;
	f_lines := strings.Split(string(f), "\n")
	re := regexp.MustCompile("^.+?\\(.+?([0-9].+?) .*?([0-9].+?) \\) PB$")
	for _, line := range f_lines {
		if strings.Contains(line, "SIGMA_BORN") {
			matches_born = re.FindAllStringSubmatch(line, -1)
		}
		if strings.Contains(line, "SIGMA_NLO") {
			matches_nlo = re.FindAllStringSubmatch(line, -1)
		}
	}
	return matches_born[0], matches_nlo[0]
}

func save_cross_sections(out string, born_vals []string, born_errs []string, nlo_vals []string, nlo_errs []string) {
	f, err := os.Create(out)
	check(err)
    defer f.Close()

	fmt.Fprintln(f, strings.Join(born_vals, `,`))
	fmt.Fprintln(f, strings.Join(born_errs, `,`))
	fmt.Fprintln(f, strings.Join(nlo_vals,	`,`))
	fmt.Fprintln(f, strings.Join(nlo_errs,	`,`))
}

func main() {
	flag.Usage = func() {
	    fmt.Fprintf(os.Stderr, "Usage: go run create_files.go 1 2 3\n")
	}   
	flag.Parse()
	parameters := flag.Args()

	base_dir := fp.Join(string(os.Getenv("HOME")), "hpair")
	hpair_in_dir := fp.Join(base_dir, "inputs/")
	hpair_out_dir := fp.Join(base_dir, "outputs/")
	create_dir(hpair_in_dir)
	create_dir(hpair_out_dir)


	// open template file
	template, err := os.ReadFile("hpair.template.in")
	check(err)

	var npars int = len(parameters)
	var sigma_born, error_born = make([]string, npars), make([]string, npars)
	var sigma_nlo, error_nlo   = make([]string, npars), make([]string, npars)

	var wg sync.WaitGroup
	for idx, tri := range parameters {
		wg.Add(1)

		go func(idx int, tri string) {
			defer wg.Done()
			
			hpair_in  := fp.Join(hpair_in_dir,  "hpair_" + tri + ".in")
			hpair_out := fp.Join(hpair_out_dir, "hpair_" + tri + ".out")
			
			create_input_file(tri, template, hpair_in)
			
			run_hpair(hpair_in, hpair_out)
			
			born, nlo := extract_cross_sections(hpair_out)

			// save data of a single goroutine
			sigma_born[idx] = born[1]
			error_born[idx] = born[2]	
			sigma_nlo[idx]  = nlo[1]
			error_nlo[idx]  = nlo[2]
		}(idx, tri)
	}
	wg.Wait()
	
	save_cross_sections(fp.Join(hpair_out_dir, "results.txt"),
		sigma_born, error_born, sigma_nlo, error_nlo)
}

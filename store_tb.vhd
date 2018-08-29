library IEEE;
use IEEE.std_logic_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

library work;
use work.common.all;

entity store_tb is
end store_tb;

architecture behavioral of store_tb is

constant clk_period : time := 10 ns;
signal clk : std_logic;
signal reset : std_logic;
signal cpuout : word;

component store is
	port (reset : in  std_logic;
				clk   : in  std_logic;
				y : out word);
end component;

begin
	 
	 u0: store port map(
				reset => reset,
				clk => clk,
				y => cpuout);
				
    proc_clock: process
    begin
        clk <= '0';
        wait for clk_period/2;
        clk <= '1';
        wait for clk_period/2;
    end process;

    proc_stimuli: process
    begin
        reset <= '1';
        wait for clk_period * 2;
        reset <= '0';
		  wait for clk_period * 16;
		  assert (cpuout = x"00000078") report "error 'cpuout' is " 
				& integer'image(to_integer(unsigned(cpuout))) severity failure;
        assert false report "success - end of simulation" severity failure;
    end process;
end architecture;

library ieee;
use ieee.std_logic_1164.all;
use ieee.std_logic_unsigned.all;

library work;
use work.common.all;

entity imem is
    port(   
        addr : in std_logic_vector(3 downto 0);
        dout : out word);
end imem;

architecture behavioral of imem is
type rom_arr is array(0 to 15) of word;

constant mem:rom_arr:=
    ( 
		x"0AA00113", -- 1301A00A              li      x2,0x0AA
		x"02100093", -- 93001002              li      x1,0x21
		x"FE20AFA3", -- A3AF20FE              sw      x2,-1(x1)
		x"022020A3", -- A3202002              sw      x2,0x21(x0)
		x"0020A0A3", -- A3A02000              sw      x2,1(x1)
		x"00000013", -- 13000000 		nop
		x"00000013", -- 13000000 		nop
		x"00000013", -- 13000000 		nop
		x"00000013", -- 13000000 		nop
		x"00000013", -- 13000000 		nop
		x"00000013", -- 13000000 		nop
		x"00000013", -- 13000000 		nop
		x"00000013", -- 13000000 		nop
		x"00000013", -- 13000000 		nop
		x"00000013", -- 13000000 		nop
		x"00000013"); -- 13000000 		nop


begin
	dout<=mem(conv_integer(addr));
end behavioral;

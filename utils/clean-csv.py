#!/usr/bin/env python3
import re
import sys


text = sys.stdin.read()
text = re.sub(r'[^\x00-\x7f]', '', text)
text = re.sub(r"([^\s'])\s+'", r"\1'", text)
text = re.sub(r"'\s+([^s'])", r"'\1", text)
text = re.sub(r"\s+/", r"/", text)
sys.stdout.write(text)

import re, sys

src = "DGSS_Model.sysml"
with open(src, "r", encoding="utf-8", newline="") as f:
    txt = f.read()
orig = txt

edits = []

# ---------- FIX① + FIX② : Initialization 状态 (删 + 标记, 统一为单个复合 entry) ----------
pat1 = (
    r'            state Initialization \{\r\n'
    r'(?:\+[^\r\n]*\r\n)+'
    r'                then action moveDropletToFDSRef;\r\n'
    r'(?:                    [^\r\n]*\r\n)+?'
    r'                then action verifyInitComplete;\r\n'
    r'(?:                    [^\r\n]*\r\n)+?'
    r'            \}'
)
repl1 = (
    '            state Initialization {\r\n'
    '                // \u2605 FIX\u2460\u2461: \u5220\u9664\u884c\u9996 + \u8865\u4e01\u6807\u8bb0; \u7edf\u4e00\u4e3a\u5355\u4e2a\u590d\u5408 entry \u52a8\u4f5c (\u4e0e Idle/Locking \u98ce\u683c\u4e00\u81f4), doc \u5185\u7f6e\u4e8e\u5404 action\r\n'
    '                entry action initSequence {\r\n'
    '                    first start;\r\n'
    '                    then action homeDGS {\r\n'
    '                        doc /* DGS \u56de\u96f6: \u6240\u6709\u5e73\u79fb/\u65cb\u8f6c\u8f74\u6267\u884c homing \u5e8f\u5217, \u786e\u8ba4\u7f16\u7801\u5668/\u9650\u4f4d\u5f00\u5173\u72b6\u6001. (\u884c\u4e3a\u5360\u4f4d) */\r\n'
    '                        assign currentState := DGSSOperatingState::initialization;     // \u2605 D5\r\n'
    '                    }\r\n'
    '                    then action moveDropletToFDSRef {\r\n'
    '                        doc /* \u63a7\u5236 DGS \u4f7f\u6db2\u6ef4\u843d\u70b9\u5bf9\u51c6 FDS \u53c2\u8003\u4f4d\u7f6e (monitorPosition), \u5b8c\u6210\u7c97\u5bf9\u51c6. */\r\n'
    '                    }\r\n'
    '                    then action verifyInitComplete {\r\n'
    '                        doc /* \u6821\u9a8c\u6240\u6709\u8f74\u4f4d\u7f6e\u5728\u671f\u671b\u533a\u95f4\u5185, \u65e0\u6545\u969c\u6807\u5fd7. */\r\n'
    '                    }\r\n'
    '                    then done;\r\n'
    '                }\r\n'
    '            }'
)
edits.append(("FIX1+2 Initialization", pat1, repl1))

# ---------- FIX③ : 属性 (uncollectableDetected 锁存 + lastDropletMissedTC + resetCommand) ----------
pat2 = (
    r'        attribute currentState : DGSSOperatingState;\r\n'
    r'(?:        //[^\r\n]*\r\n)+'
    r'        attribute uncollectableDetected : Boolean;\r\n'
    r'(?:        //[^\r\n]*\r\n)+'
    r'        attribute faultFlag    : Boolean = uncollectableDetected;'
)
repl2 = (
    '        attribute currentState : DGSSOperatingState;\r\n'
    '        // \u2605 FIX\u2462 (\u95ed\u5408 TODO(D7)): \u4e0d\u53ef\u6536\u96c6\u68c0\u6d4b\u9501\u5b58. \u521d\u503c false; \u8fd0\u884c\u6001\u7531 do \u52a8\u4f5c\u7f6e\u4f4d (\u89c1 monitorStability/trackReference),\r\n'
    '        //   Fault \u590d\u4f4d\u65f6\u7531 faultToInit \u7684 effect \u6e05\u96f6. \u7f6e\u4f4d\u8bed\u4e49: pointInsideTCAperture(...)==false \u2228 TC-miss\u2265N \u2228\r\n'
    '        //   FDS-loss \u2228 \u8d85\u884c\u7a0b, \u8fd9\u4e9b\u6761\u4ef6\u7ecf FDS/\u51e0\u4f55\u5224\u5b9a\u6c47\u5165 lastDropletMissedTC.\r\n'
    '        attribute uncollectableDetected : Boolean = false;\r\n'
    '        // \u2605 FIX\u2462: FDS/\u51e0\u4f55\u5224\u5b9a\u8f93\u51fa \u2014 \u6700\u8fd1\u4e00\u9897\u6db2\u6ef4 TC \u7a7f\u8d8a\u70b9\u662f\u5426\u843d\u5728\u5f00\u53e3\u5916 (\u5916\u90e8\u68c0\u6d4b\u8f93\u5165\u4fe1\u53f7).\r\n'
    '        attribute lastDropletMissedTC : Boolean = false;\r\n'
    '        // \u2605 FIX\u2462: \u64cd\u4f5c\u5458\u624b\u52a8\u590d\u4f4d\u6307\u4ee4 (\u5916\u90e8\u8f93\u5165; \u89c1 Safety-1 \u5e94\u5bf9\u63aa\u65bd\u7b2c 7 \u6b65). \u89e6\u53d1 Fault\u2192Initialization \u5e76\u6e05\u9501\u5b58.\r\n'
    '        attribute resetCommand : Boolean = false;\r\n'
    '        // \u2605 D4: faultFlag \u7ed1\u5b9a\u5230\u68c0\u6d4b\u9501\u5b58 (\u539f `= false` \u6052\u5047 \u2192 \u6545\u969c\u8fc1\u79fb\u6c38\u4e0d\u53ef\u8fbe).\r\n'
    '        attribute faultFlag    : Boolean = uncollectableDetected;'
)
edits.append(("FIX3 attributes", pat2, repl2))

# ---------- FIX③ : Idle do-action 折叠 (置位 producer) ----------
pat3 = (
    r'                do action monitorStability;\r\n'
    r'(?:                    [^\r\n]*\r\n)+'
)
repl3 = (
    '                do action monitorStability {\r\n'
    '                    doc /* \u6301\u7eed\u4ece FDS \u8bfb\u53d6\u6db2\u6ef4\u7a33\u5b9a\u6027\u6d4b\u91cf\u7ed3\u679c, \u4f46\u4e0d\u95ed\u73af\u4fee\u6b63. */\r\n'
    '                    assign uncollectableDetected := uncollectableDetected or lastDropletMissedTC;   // \u2605 FIX\u2462: \u9501\u5b58\u7f6e\u4f4d\r\n'
    '                }\r\n'
)
edits.append(("FIX3 Idle.do", pat3, repl3))

# ---------- FIX③ : ClosedLoopLocking do-action 折叠 (置位 producer) ----------
pat4 = (
    r'                do action trackReference;\r\n'
    r'(?:                    [^\r\n]*\r\n)+'
)
repl4 = (
    '                do action trackReference {\r\n'
    '                    doc /* \u5b9e\u65f6\u8ba1\u7b97\u6db2\u6ef4\u4f4d\u7f6e\u8bef\u5dee, \u8f93\u51fa DGS \u4fee\u6b63\u6307\u4ee4. */\r\n'
    '                    assign uncollectableDetected := uncollectableDetected or lastDropletMissedTC;   // \u2605 FIX\u2462: \u9501\u5b58\u7f6e\u4f4d\r\n'
    '                }\r\n'
)
edits.append(("FIX3 Locking.do", pat4, repl4))

# ---------- FIX③ : faultToInit 复位 (清锁存 effect) ----------
pat5 = (
    r'            transition faultToInit\r\n'
    r'                first Fault\r\n'
    r'                accept after \(1 \[s\]\)[^\r\n]*\r\n'
    r'                if \(not faultFlag\)\r\n'
    r'                then Initialization;'
)
repl5 = (
    '            transition faultToInit                                              // \u2605 FIX\u2462: \u590d\u4f4d\u7531\u64cd\u4f5c\u5458\u6307\u4ee4\u89e6\u53d1, effect \u6e05\u9501\u5b58\r\n'
    '                first Fault\r\n'
    '                if (resetCommand)\r\n'
    '                do action clearFaultLatch {\r\n'
    '                    assign uncollectableDetected := false;     // \u6e05\u9664\u6545\u969c\u9501\u5b58\r\n'
    '                    assign resetCommand          := false;     // \u590d\u4f4d\u6307\u4ee4\u81ea\u6e05\u9664 (\u5355\u6b21\u89e6\u53d1)\r\n'
    '                }\r\n'
    '                then Initialization;'
)
edits.append(("FIX3 faultToInit", pat5, repl5))

# ---- apply with uniqueness checks ----
for name, pat, repl in edits:
    n = len(re.findall(pat, txt))
    if n != 1:
        print(f"[ABORT] {name}: matched {n} times (expected 1)"); sys.exit(1)
    txt = re.sub(pat, lambda m: repl, txt, count=1)
    print(f"[OK]    {name}: applied")

with open(src, "w", encoding="utf-8", newline="") as f:
    f.write(txt)

print("\n--- post-checks ---")
print("leftover '+'-prefixed code lines:", len(re.findall(r'\r\n\+[^\r\n]*', txt)))
for tok in ["lastDropletMissedTC", "resetCommand", "clearFaultLatch", "initSequence", "FIX\u2460\u2461", "FIX\u2462"]:
    print(f"  contains {tok!r}:", tok in txt)
print("line-count orig vs new:", orig.count('\n'), txt.count('\n'))

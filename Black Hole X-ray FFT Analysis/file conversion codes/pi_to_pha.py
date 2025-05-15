from astropy.io import fits

# 読み取り
hdul = fits.open("ad15606000s010102_1.pi")

# PHAフォーマットの条件確認
for hdu in hdul:
    if hdu.header.get("EXTNAME") == "SPECTRUM":
        print("PHA互換の拡張子が存在します。")

# ファイル名変更のみなら
hdul.writeto("ad15606000s010102_1.pha", overwrite=True)

# 完了
print("piファイルをphaファイルとして書き出しました。")

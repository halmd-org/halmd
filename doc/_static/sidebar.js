$(window).scroll(function() {
    var margin = 10;
    var content_top = $(".content-wrapper").position().top;
    var header_height = $(".logo").outerHeight() + $(".header-wrapper").outerHeight() + margin;
    var content_height = $(".document").outerHeight();
    var sidebar_height = $(".sidebar").outerHeight();
    var scroll_top = $(window).scrollTop();


    // fix the sidebar to the header if it is visible or the content too short
    if (scroll_top < header_height || content_height < sidebar_height) {
        $(".sidebar").css({position: 'static', top: Math.max(margin, -content_top), 'margin-left' :0 });
    }
    else {
        // make sure we dont move the sidebar over the footer
        var footer_break = $(document).height() - (sidebar_height + $(".footer-wrapper").outerHeight() + margin);
        var top = footer_break - scroll_top;
        if (top > margin) top = margin;
        $(".sidebar").css({position: 'fixed', top : top, 'margin-left' : $(".content").width()-$(".sidebar").width() });
    }
});
$(window).resize = $(window).scroll;